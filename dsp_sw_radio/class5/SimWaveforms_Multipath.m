num_symbols = 2400;
qpsk = 4;
qam16 = 16;
oversamp = 4;
alpha = 0.25;
sym_rate = 1e6;
sym_length = 46;
[tx_shaped,QAM_Const] = rrc_qam(num_symbols, qpsk, oversamp, alpha, sym_length);
figure()
plot(real(tx_shaped),imag(tx_shaped))

eyediagram(real(tx_shaped),2*oversamp, 2);
title("Eye Diagram of RRC Tx Waveform (Real)");

coeff = psCoeff(oversamp, alpha, 'rrc', sym_length);

rx_shaped = filter(coeff, 1, tx_shaped);

transient = (sym_length/2 + 1) * oversamp;

eyediagram(rx_shaped(transient:end),2*oversamp, 2);
title("Eye Diagram of Rx Waveform - After 2nd RRC Filter");


% Initial EVM   (Note decisions below are QPSK)
offset = 3;  % offset to select correct closest decision locations out of the 4 samples/symbol
             % (adjust between 0 and 3 for minimum EVM)

decisionI = (round(real(rx_shaped(transient + offset : 4 : end))/2 + 0.5) - 0.5)/2;
decisionQ = (round(imag(rx_shaped(transient + offset : 4 : end))/2 + 0.5) - 0.5)/2;
decisions = decisionI + 1i * decisionQ;

evm_value = evm(rx_shaped(transient + offset : 4 : end), decisions);
evm_dB = db(evm_value);

fprintf('EVM at transmitter = %.1f dB\n', evm_dB);

% Confirm scaling is not limiting EVM
mean_abs_decisions = mean(abs(decisions));
mean_abs_rx_shaped = mean(abs(rx_shaped(transient + offset : 4 : end)));

fprintf('Mean of absolute decisions: %.4f\n', mean_abs_decisions);
fprintf('Mean of absolute rx_shaped samples: %.4f\n', mean_abs_rx_shaped);

scaling_error_evm = db(mean_abs_decisions - mean_abs_rx_shaped);
fprintf('EVM as limited by scaling: %.1f dB\n', scaling_error_evm);
% this is the EVM estimate contribution due to scaling error alone

rx_scaled = rx_shaped * mean_abs_decisions / mean_abs_rx_shaped;
scaled_error_evm = db(evm(rx_scaled(transient + offset : 4 : end), decisions));
fprintf('EVM after correcting scale error: %.1f dB\n', scaled_error_evm);

scaled_error_limit = db(mean_abs_decisions - mean(abs(rx_scaled(transient + offset : 4 : end))));
fprintf('EVM as limited by scaling after correcting scale error: %.1f dB\n', scaled_error_limit);


%% TX Spectrum
fclk = sym_rate*oversamp;
figure()
% using cusomized fft module imported earlier
[Pxx,f] = win_fft(tx_shaped/max(tx_shaped), sym_rate*oversamp);
plot(f-fclk/2, 20*log10(abs(fftshift(Pxx))), 'LineWidth', 2);
title('NCO Output Spectrum')
title('QPSK Tx Spectrum');

%% Multipath Distortion Demonstration
% model a multipath channel with five dominant but similarly distributed paths
cir = [0.019*exp(1i*0.2), 0.021*exp(-1i*0.24), 0.08*exp(1i*2.4), 0.032*exp(-1i*0.14), -0.032*exp(-1i*0.14)];

[H, w] = freqz(cir, 1, 'whole');
w = w - pi;

figure;
subplot(2,1,1);
plot(w, db(fftshift(H)));
grid on;
title('Channel Frequency Response');
ylabel('Magnitude [dB]');

subplot(2,1,2);
plot(w, unwrap(fftshift(angle(H))) * 360 / (2 * pi));
grid on;
ylabel('Phase [deg]');
xlabel('Normalized Radian Frequency');

% Pass tx waveform through channel
rx_distort = filter(cir, 1, tx_shaped);

[H, w] = freqz(cir, 1, 'whole');
w = w - pi;

%% channel group delay
[gd1, w_gd] = grpdelay(cir, 1, 512, 'whole');
w_gd = w_gd - pi;

figure;
plot(w_gd, fftshift(gd1));
grid on;
title('Channel Group Delay');
ylabel('Delay [samples]');
xlabel('Normalized Radian Frequency');
axis([-pi, pi, 0, 4]);

%% Eye Diagram after 2nd RRC Filter - Unequalized Multipath Rx
rx_pre_equal = filter(coeff, 1, rx_distort);

% Adjusting the transient for the increased oversampling factor (64)
transient_equalized = transient * 64 / oversamp;

% Generating an eye diagram for the distorted waveform after the second RRC filter
eyediagram(rx_pre_equal(transient_equalized:end), 2*oversamp, 2);

title("Eye Diagram of Distorted Waveform after 2nd RRC Filter");

%% Rx Spectrum (Unequalized)
figure()
% using cusomized fft module imported earlier
[Pxx_rx_distort,F_rx_distort] = win_fft(rx_distort, sym_rate*oversamp);
plot(F_rx_distort-fclk/2, 20*log10(abs(fftshift(Pxx_rx_distort))), 'LineWidth', 2);
title('QPSK Rx Spectrum Prior to Equalization');
xlabel('Normalized Radian Frequency');

%% 
figure;
plot(real(rx_distort(1:500)) / std(rx_distort), 'DisplayName', 'Rx');
hold on;
plot(real(tx_shaped(1:500)) / std(tx_shaped), 'DisplayName', 'Tx');
hold off;

title('Tx and Rx Waveform');
xlabel('Samples');
legend;

% Omitting the first 100 samples of the sequence, which is zero due to filter loading delays
rx = rx_distort(101:end);
tx = tx_shaped(101:end);

%% LMS Equalizer 
omit = 100;    % initial samples to exclude
shift = 0;     % number of samples to shift dominant equalizer tap to the right

ntaps = 75;    % increase to see the number of significant coefficients needed below
depth = ntaps * 80;
delay = floor(ntaps / 2);
A = convmtx(rx_distort(omit + shift : omit + shift + depth), ntaps);
R = conj(A.') * A;
X = [zeros(1, delay), tx_shaped(omit : omit + depth), zeros(1, ceil(ntaps / 2) - 1)];
ro = conj(A.') * X.';

equal_coeff = R \ ro;


%% Plotting magnitudes of equalizer coefficients to show dominant taps
% If the dominant taps favor the left or right side, change the variable shift accordingly
% (increase shift positive to shift dominant tap to the right...)
figure;
stem(abs(equal_coeff));
title('Magnitude of Equalizer Coefficients');
ylabel('Magnitude');

%% Equalizer frequency response
[H_equalizer, w_equalizer] = freqz(equal_coeff, 1, 'whole');
w_equalizer = w_equalizer - pi;

figure;
subplot(2,1,1);
plot(w_equalizer, db(fftshift(H_equalizer)));
grid on;
title('Equalizer Frequency Response');
ylabel('Magnitude [dB]');

subplot(2,1,2);
plot(w_equalizer, unwrap(fftshift(angle(H_equalizer))) * 360 / (2 * pi));
grid on;
ylabel('Phase [deg]');
xlabel('Normalized Radian Frequency');

%% Equalizer Group Delay
[ gd_equalizer, w_equalizer] = grpdelay(equal_coeff, 1, 512,'whole');
w_equalizer = w_equalizer - pi;

figure;
plot(w_equalizer, fftshift(gd_equalizer));
grid on;
title('Equalizer Group Delay');
ylabel('Delay [samples]');
xlabel('Normalized Radian Frequency');
%% corrected group delay
figure;
plot(w_equalizer, fftshift(gd1 + gd_equalizer));
grid on;
title('Corrected Group Delay');
ylabel('Delay [samples]');
xlabel('Normalized Radian Frequency');
% axis([-pi, pi, ntaps/2-5, ntaps/2+5]);

%% pass distorted waveform through equalizer
rx_recovered = filter(equal_coeff,1,rx_distort);

% Adjusting the transient for the increased oversampling factor (64)
transient_equalized = transient * 64 / oversamp;

% Generating an eye diagram for the distorted waveform after the second RRC filter
eyediagram(rx_recovered(transient_equalized:end), 2*oversamp, 2);

title("Eye Diagram of Equalized Waveform");

%% Eye Diagram after 2nd RRC Filter - Equalized Multipath Rx
rx_equalized = filter(coeff, 1, rx_recovered);

% Adjusting the transient for the increased oversampling factor (64)
transient_equalized = transient * 64 / oversamp;


eyediagram(rx_equalized, 2 * oversamp, 2);
title('Eye Diagram of Equalized Waveform');

%%  EVM after equalization
% (Typical EVM requirements listed below (actual requirement depends on coding used and BER requirements):
% QPSK: -10 to -13 dB
% 16-QAM: -16 to -19 dB
% 64-QAM: -22 to -27 dB

offset = 4;  % offset to select correct closest decision locations out of the 4 samples/symbol
             % set from eye diagram above, choosing sample number at the decision location

decisionI = (round(real(rx_equalized(transient + offset : 4 : end))/2 + 0.5) - 0.5)/2;
decisionQ = (round(imag(rx_equalized(transient + offset : 4 : end))/2 + 0.5) - 0.5)/2;
decisions = decisionI + 1j * decisionQ;

evm_value_equalized = evm(rx_equalized(transient + offset : 4 : end), decisions);
evm_dB_equalized = db(evm_value_equalized);

fprintf('EVM after equalization = %.1f dB\n', evm_dB_equalized);

%% Equalized spectrum 
figure()
% using cusomized fft module imported earlier
[Pxx_rx_equalized,F_rx_equalized] = win_fft(rx_equalized/max(rx_equalized), sym_rate*oversamp);
plot(F_rx_equalized-fclk/2, 20*log10(abs(fftshift(Pxx_rx_equalized))), 'LineWidth', 2);
title('QPSK Rx Spectrum After Equalization');
xlabel('Normalized Radian Frequency');


%% By swapping input and output we can estimate the channel rather than equalize for it:
% Compute channel coefficients of estimate:
omit = 100;    % Initial samples to exclude
shift = 0;     % Number of samples to shift dominant equalizer tap to the right

ntaps = 25;    % Increase to see the number of significant coefficients needed below
depth = ntaps * 80;
delay = floor(ntaps / 2);
A = convmtx(tx_shaped(omit + shift : omit + shift + depth), ntaps);
R = conj(A.') * A;
X = [zeros(1, delay), rx_distort(omit : omit + depth), zeros(1, ceil(ntaps / 2) - 1)];
ro = conj(A.') * X.';

channel_coeff = R \ ro;

figure()
plot(abs(channel_coeff))


%%
[H_actual, w_actual] = freqz(cir, 1, 'whole');
[H_estimated, w_estimated] = freqz(channel_coeff, 1, 'whole');
w_actual = w_actual - pi;
w_estimated = w_estimated - pi;

figure;
plot(w_actual, 20*log10(abs(fftshift(H_actual))));
hold on;
plot(w_estimated, 20*log10(abs(fftshift(H_estimated))));
hold off;

xlabel('Normalized Radian Frequency');
ylabel('dB');
title('Channel Estimation Using Wiener-Hopf Equations');
legend("Actual Channel", "Estimated Channel");
grid on;
%% 
X = fft(tx_shaped(omit : omit + depth), 2 * depth);
Y = fft(rx_distort(omit : omit + depth), 2 * depth);

figure;
plot(20 * log10(abs(X)));
hold on;
plot(20 * log10(abs(Y)));
hold off;

figure;
plot(20 * log10(abs(Y ./ X)));
%% 
channel = ifft(Y ./ X);
figure;
plot(abs(channel));

figure;
[h_estimated,w_estimated] = freqz(channel, 1, 'whole');
w_estimated = w_estimated - pi;
plot(w_estimated, 20*log10(abs(fftshift(h_estimated))));
hold on;
plot(w_actual, 20*log10(abs(fftshift(H_actual))));
hold off;

xlabel('Normalized Radian Frequency');
ylabel('dB');
title('Channel Estimation Using IFFT(FFT(y)/FFT(x))');
legend("Estimated Channel", "Actual Channel");
grid on;





%% repeating above equalization with a carrier offset condition
% given equalizer span of 5 symbols with 4 samples per symbol T= 20
% a carrier offset in normalized to the symbol rate
fracf = 0.0001;   % normalized frequency offset of carrier relative to symbol rate 
                  %          offset   EVM
                  %          0.0001  -34.4 dB (100 Hz offset for 1MHz symbol rate)
                  %          0.0002  -30.8 dB (200 Hz)
                  %          0.0003  -27.3 dB (300 Hz)
                  %          0.0004  -24.2 dB (400 Hz)
                  %          0.0005  -13.8 dB (500 Hz) 
            
n = 0:length(rx_pre_equal)-1;

rx_offset1 = rx_distort .* exp(1i * 2 * pi * n / oversamp * fracf);

% Compute equalizer coefficients
A = convmtx(rx_offset1(omit + shift : omit + shift + depth), ntaps);
R = A' * A;
X = [zeros(1, delay), tx_shaped(omit : omit + depth), zeros(1, ceil(ntaps / 2) - 1)];
ro = A' * X.';

equal_coeff = R \ ro;

% Pass distorted waveform through equalizer
rx_recovered = filter(equal_coeff, 1, rx_offset1);

% Remove carrier offset
rx_after_car_recovery = rx_recovered .* exp(-1i * (2 * pi * n / oversamp * fracf));

% 2nd RRC filter and phase correction
rx_equalized = filter(coeff, 1, rx_after_car_recovery);

% Correct any residual phase offset with decision-directed phase detector
decisionI = (round(real(rx_equalized) / 2 + 0.5) - 0.5)/2;
decisionQ = (round(imag(rx_equalized) / 2 + 0.5) - 0.5)/2;
decisions = decisionI + 1i * decisionQ;

offset = 2;    % Decision location (confirm matches eye diagram or update)
start = 4 * 27;  % Eliminate initial transients

ph_errors = angle(decisions(offset:4:end) .* conj(rx_equalized(offset:4:end)));
disp(ph_errors);
disp(['Residual phase error = ' num2str(mean(ph_errors))]);

% rx_after_phase_correction =  rx_after_car_recovery;
rx_after_phase_correction = rx_equalized .* exp(1i * mean(ph_errors(start:end)));

figure;
plot(ph_errors(start:end));
title('Residual Constellation Rotation (Corrected)');
xlabel('Samples');
ylabel('Radians');
% axis([0, 2500, -pi, pi]);


eyediagram(rx_after_phase_correction(transient*64/oversamp:end), 2 * oversamp, 2);
title('Eye Diagram of Equalized Waveform');

%%
offset = 1;

figure;
plot(real(rx_after_phase_correction), imag(rx_after_phase_correction));
hold on;
plot(real(rx_after_phase_correction(start+offset:4:end)), imag(rx_after_phase_correction(start+offset:4:end)), 'r.');
hold on
plot(decisionI, decisionQ, '+');
hold off;

