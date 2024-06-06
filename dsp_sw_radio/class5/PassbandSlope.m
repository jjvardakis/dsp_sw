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

% create impulse response for bandlimited passband slope as our "channel filter"
ntaps = 99;
sym_length = 46;
oversamp = 4;
% bandlimited in frequency is a Sinc in time, and we'll window this with a Kaiser window to reduce truncation error:
fc = 2/(2*oversamp);  % one-sided bandwidth as a fraction of the sampling rate, 1/4 would be a "half-band" filter    
n = -(ntaps-1)/2:(ntaps-1)/2;
coeff_real = 2 * fc * sinc(2*fc*n);  % normalized channel gain = 1
disp(coeff_real);

% The complex impulse response will result in a passband slope
slope = -0.2;  % slope is the change in magnitude versus the change in frac frequency in rad/sample (2pi = sampling rate)
coeff_imag = 1j * slope * 2 * pi / (100 * fc) * dsinc(2*fc*n);

% Another approach to get passband slope:
% I decomposed freq response into a rect + integral of rect to get passband slope
% - The inverse FT of rect is a real Sinc in time
% - The inverse FT of the integral of the rect is Sinc/(2pi t) given integral of X(f)df ---> j x(t)/(2pi f)
% With this approach, the slope will continue in the stopband since the slope portion of the frequency response 
% would be actually be an integration of the rect, multiplied by the rect. Since there is no waveform
% energy in the stop band, this doesn't matter but if we did care we would 
% have to convolve the simpler time domain impulse response result j x(t)/(2pi f) with another Sinc. 
% I did not do that in the code commented out below, so if we used it we'll see an additional 
% slope in the stopband using this approach.

% num = slope * pi / fc * 1j * coeff_real;  % scales slope to magnitude vs rad/sample
% den = 2 * pi * n;
% dividing without divide by zero error: 
% coeff_imag = num ./ den;
% coeff_imag(isnan(coeff_imag) | isinf(coeff_imag)) = 0;

channel = coeff_real + coeff_imag;
% window result
channel_win = channel .* kaiser(length(channel), 6).';

figure;
subplot(2,1,1);
stem(n, coeff_real, 'MarkerSize', 3);
title('Bandlimited Channel Impulse Response - Real');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
stem(n, imag(coeff_imag), 'MarkerSize', 3);
title('Bandlimited Channel Impulse Response - Imag');
xlabel('Time (samples)');
ylabel('Amplitude');
grid on;


%% interpolated impulse response to see ISI impact directly:
m = linspace(-(ntaps-1)/2, (ntaps-1)/2 + 1, 199);

real_sinc = 2 * fc * sinc(2 * fc * m);
imag_dsinc = 1j * slope / 100 * oversamp * 2 * pi * dsinc(2 * fc * m);

figure;
plot(real_sinc);
hold on;
plot(imag(imag_dsinc));
hold off;
grid on;

% The Sinc impulse response is a Nyquist pulse that follows the Nyquist ISI criterion and results in going through
% zero at integer multiples of the symbol rate.
% We see in the plot below and earlier descriptions above that in order to get the passband slope, we needed
% a derivative of that sinc as the imaginary component of the modified impulse response. The derivative has
% zeros at the peaks of the Sinc and therefore will not be zero at the symbol boundary locations: ISI!

%% Channel frequency response
[h, w] = freqz(channel_win,1,'whole');

figure;
subplot(2,1,1);
plot(w - pi, abs(fftshift(h)));
xlabel('Frequency');
ylabel('Magnitude');
grid on;

subplot(2,1,2);
plot(w - pi, angle(fftshift(h)));
xlabel('Frequency');
ylabel('Phase');
grid on;

%% Showing effect on constellation
rx_distortion = filter(channel, 1, rx_shaped);
sym_length = 46;
transient = (sym_length/2 + 1) * oversamp;

eyediagram(rx_distortion(transient*64/oversamp+1:end), 2 * oversamp, 2);
title("Eye Diagram of Waveform after Channel Distortion");

% And here we see the ISI mentioned earlier. The real portion of the zero-ISI symbols creates
% ISI in the imaginary, and the imaginary portion of the zero-ISI symbols creates ISI in the real.
% So in additional detail to ISI, it is an ISI leakage from real to imaginary. (So then a BPSK waveform
% alone should maintain zero-ISI)


%% Constellation after distortion
start_offset = 50;
timing_offset = 2;    % enter sample offset from eye diagram above
i_after = real(rx_distortion(start_offset * oversamp + timing_offset:4:end));
q_after = imag(rx_distortion(start_offset * oversamp + timing_offset:4:end));

disp(i_after(1:10));

% Constellation before channel distortion:
start_offset = 25;
timing_offset = 5;
i_before = real(rx_shaped(start_offset * oversamp + timing_offset:4:end));
q_before = imag(rx_shaped(start_offset * oversamp + timing_offset:4:end));

figure;
plot(i_after, q_after, 'g.', 'DisplayName', 'After Distortion');
hold on;
plot(i_before, q_before, 'r.', 'DisplayName', 'Before Distortion');
hold off;
legend;
grid on;
axis equal;

%% Actual EVM measured:
decisionI = (round(i_after/2+0.5)-0.5)/2;
decisionQ = (round(q_after/2+0.5)-0.5)/2;
decisions = decisionI + 1j * decisionQ;
received = i_after + 1j * q_after;

EVM = 20 * log10(std(received - decisions) / std(decisions));
fprintf('EVM = %.2f dB\n', EVM);

%% BPSK case: no distortion!
% All the distortion from passband slope is in the imaginary axis so does not show up in a BPSK signal 
% For ISI to occur there would also need to be a quadrature signal to create a real distortion.

rx_distortion_bpsk = filter(channel, 1, real(rx_shaped));

eyediagram(rx_distortion_bpsk(transient*64/oversamp+1:end), 2 * oversamp, 2);
title("Eye Diagram of Waveform after Channel Distortion");


%% Create random QPSK waveform

% Set the seed for reproducibility
rng(42);

% Create random data source
num_symb = 10000;

% Undistorted Tx waveform with no quadrature error as raised cosine pulse shaping
% (instead of RRC) mimicking what would occur through the cascade of Tx and Rx

% Assuming rc_qam function definition (replace with the actual implementation)
tx_no_error = rc_qam(num_symb, qpsk, oversamp, alpha, sym_length);

samps = plotConstellation(tx_no_error, 42, 'QPSK Waveform with no IQ Erro', 4, 10);

%% Induce IQ error in amplitude and phase and demonstrate adaptive correction technique

% Consider an IQ mixer such as https://www.analog.com/media/en/technical-documentation/data-sheets/5589f.pdf
% that contains an analog quadrature splitter on the local oscillator,
% the output of the splitter provides a sine and cosine reference of the L.O.
% If the two are not in perfect quadrature, we can consider the cosine (in phase) output
% to be the "reference phase" and the sine (quadrature) output to be in error such
% that some of the Q modulation will leak into the I mixer output.

phase_err = 0.2;          % introduced quadrature error (typically <0.2 radians)
rot = 0.6;                % introduced arbitrary rotation (change to different values)

% Distorted waveform with quadrature error and rotation (Phase Offset)
tx_phase_error = (real(tx_no_error) + sin(phase_err) * imag(tx_no_error) + 1j * imag(tx_no_error)) * exp(1j * rot);

samps = plotConstellation(tx_phase_error, 42, 'QPSK Waveform with Quadrature Error and Phase Offset', 4, 10);

%% Plot samples and decisions showing approach to remove Phase Offset

% Normalize constellation 
sampsNorm = samps / std(samps);

decisions = (2 * (real(sampsNorm) > 0) - 1) / sqrt(2) + 1i * (2 * (imag(sampsNorm) > 0) - 1) / sqrt(2);

figure;
plot(real(sampsNorm), imag(sampsNorm), 'r.');
axis equal;
ax = axis;
grid on;

hold on;
plot(real(decisions), imag(decisions), 'b+');
axis([-2, 2, -2, 2]);
title('Actual Samples and Decisions');
xlabel('Real');
ylabel('Imag');
hold off;

% Alternate phase estimator using quadrature samples and sign of real:
rotCorrect = sign(real(sampsNorm)) .* (imag(sampsNorm) - imag(decisions));

figure;
plot(rotCorrect);
title('Measured Phase Error using $sign(I)(Q-\hat Q)$', 'Interpreter', 'latex');
xlabel('Sample Number');
ylabel('Error');

%% Correction loop to remove rotation
% (Type 1 first-order loop, okay for removing a static phase offset, but to 
% remove frequency offset, a Type 2 Loop would be a better choice)

accumQ = 0;
k = 0.5;

% Save results of loop simulation for plotting
errs = [];
derotatedSamples = [];

for idx = 1:numel(sampsNorm)
    sample = sampsNorm(idx);
    sample_rotated = sample * exp(-1i * accumQ);
    derotatedSamples(end + 1) = sample_rotated;
    decision = (2 * (real(sample_rotated) > 0) - 1) / sqrt(2) + 1i * (2 * (imag(sample_rotated) > 0) - 1) / sqrt(2);
    
    % New phase detector that will drive Q difference to 0 as desired for subsequent quadrature removal
    errQ = sign(real(sample_rotated)) * (imag(sample_rotated) - imag(decision));
    
    % Phase detector given by ^I*Q - I*^Q will leave a rotation offset due to quadrature imbalance
    % Although this may still be better for higher-order QAM due to consistent error result for decisions closer to 
    % imaginary axis, and when combined with quadrature correction in loop should still converge. 
    % errQ = real(decision) * imag(sample_rotated) - imag(decision) * real(sample_rotated);
    
    errs(end + 1) = errQ;
    accumQ = accumQ + k * errQ;
end

derotatedSamples = derotatedSamples.';
figure;
plot(real(errs));
title('Measured Rotation Error vs Sample in Loop');
xlabel('Sample');
ylabel('Error');

figure;
plot(real(derotatedSamples), imag(derotatedSamples), 'r.');
axis equal;
grid on;
title('Constellation after Rotation Removal - Only Quadrature Error Remains');
xlabel('Real');
ylabel('Imag');

%% Correlate I to Q to determine phase error
N = length(derotatedSamples);
corr = dot(imag(derotatedSamples), real(derotatedSamples)) / (N * std(imag(derotatedSamples)) * std(real(derotatedSamples)));
disp(['Correlation between I and Q: ', num2str(corr)]);

%% Note: The above computation shows the underlying operations but can be computed
% directly using corrcoef in MATLAB. In practice, a typical approach would simply use
% the dot product when a result proportionate to the phase error can be used rather
% than an accurate measurement (such as in a correction loop)

corr_matrix = corrcoef(imag(derotatedSamples), real(derotatedSamples))

%% Below shows the result for the waveform with no quadrature error

corr_matrix_no_error = corrcoef(imag(tx_no_error), real(tx_no_error))

%% Using the measured phase error to correct the constellation

tx_corrected_Q = imag(derotatedSamples);
tx_corrected_I = real(derotatedSamples) - imag(derotatedSamples) * corr;

tx_corrected = tx_corrected_I + 1i * tx_corrected_Q;

%% Plot the constellation after rotation and quadrature correction
figure;
plot(tx_corrected_I, tx_corrected_Q, 'r.');
axis equal;
grid on;
title('Constellation after Rotation And Quadrature Correction');
xlabel('Real');
ylabel('Imag');


%%  Quadrature Correction Loop Demo
beta = 0;
k = 0.15;    % high loop bandwidth for acquisition
alpha = 0.97;
k2 = 0.0003;

acqsamps = 50;

nsamps = floor(log(k2/k) / log(alpha));
Iout = [];
Qout = [];
betas = [];
errs = [];

disp(['Number of samples in gear shift transition = ', num2str(nsamps)]);

for count = 1:length(derotatedSamples)

    % Blind gear shifting for fast acquisition
    if count > acqsamps
        k = k * alpha;   % transition samples. Will transition in log(0.00015/0.025)/log(0.95) samples
    end
    if count > acqsamps + nsamps
        k = k2;   % narrow loop bandwidth for tracking
    end
    
    % Apply correction
    sample = derotatedSamples(count);
    Q = imag(sample);
    I = real(sample) - beta * imag(sample);
    Iout = [Iout; I];
    Qout = [Qout; Q];
    err = I * Q;
    errs = [errs; err];
    beta = beta + k * err;
    betas = [betas; beta];
end

figure()
plot(betas)
title('Convergence of Beta in Loop')
xlabel('Sample Number')
ylabel('Beta')
grid on

% Plot the constellation after quadrature correction loop
figure;
plot(Iout(4001:end), Qout(4001:end), 'r.'); % Showing samples fully converged
axis equal;
grid on;
title('Constellation After Quadrature Correction Loop');

%% alternative implementation - Quadrature Correction Loop Demo
beta = 0;
k = 0.2;
k2 = 0.002;
alpha = 0.96;

acqsamps = 20;

nsamps = floor(log(k2/k) / log(alpha));
disp(['Number of samples in gear shift transition = ', num2str(nsamps)]);

Iout = [];
Qout = [];
betas = [];
errs = [];

for count = 1:length(derotatedSamples)
    % Blind gear shifting for fast acquisition
    if count > acqsamps
        k = k * alpha;  % Transition samples. Will transition in ln(k2/k)/ln(alpha) samples
    end
    if count > acqsamps + nsamps
        k = k2;  % Narrow loop bandwidth for tracking
    end
    
    % Apply quadrature correction
    sample = derotatedSamples(count);
    Q = imag(sample);
    I = real(sample) - beta * imag(sample);
    Iout = [Iout; I];
    Qout = [Qout; Q];
    err = sign(Q) * (I - sign(I) / sqrt(2));  % This is sign(Q)(I-^I)
    errs = [errs; err];
    beta = beta + k * err;
    betas = [betas; beta];
end


figure()
plot(betas)
title('Convergence of Beta in Loop')
xlabel('Sample Number')
ylabel('Beta')
grid on

% Plot the constellation after quadrature correction loop
figure;
plot(Iout(4001:end), Qout(4001:end), 'r.'); % Showing samples fully converged
axis equal;
grid on;
title('Constellation After Quadrature Correction Loop');
%% 
snr = 14;       % SNR in dB (reduce to test performance in lower SNR conditions)

nr = 10^(-snr/10)/sqrt(2);     % Noise ratio on I and Q

noise = nr * (randn(size(sampsNorm)) + 1i * randn(size(sampsNorm)));

alpha = 0.96;
beta = 0;
theta = 0;
k_acq = 0.25;
k_track = 0.005;

acqsamps = 30;

nsamps = floor(log(k_track/k_acq) / log(alpha));
disp(['Number of samples in gear shift transition = ', num2str(nsamps)]);

Iout = [];
Qout = [];
betas = [];
errsRot = [];
errsQuad = [];
thetas = [];

acqDetectBuffer = [];
acqState = 1;
k = k_acq;
acqCounter = acqsamps;

for count = 1:length(sampsNorm)
    % Acquisition detect for gear shift
    if acqState == 1    % Acquisition
        if acqCounter < 0
            if mean(acqDetectBuffer) < 4 * std(acqDetectBuffer)
                acqState = 0;
                k = k_track;
                disp(['Tracking mode at sample # ', num2str(count), '!']);
            end
        else
            acqCounter = acqCounter - 1;
        end
    end
    if acqState == 0    % Tracking
        if mean(acqDetectBuffer) > 4 * std(acqDetectBuffer)
            acqState = 1;
            k = k_acq;
            acqCounter = acqsamps;
            disp(['Back to acquisition mode at sample # ', num2str(count), '!']);
        end
    end
    sample = sampsNorm(count) + noise(count);
    % Apply corrections
    sample_rotated = sample .* exp(-1i*theta);
    Q = imag(sample_rotated);
    I = real(sample_rotated) - beta * imag(sample_rotated);

    decision = sign(I)/sqrt(2) + 1i * sign(Q)/sqrt(2);

    errRot = sign(real(sample_rotated)) * (imag(sample_rotated) - imag(decision));  % Drive Q delta to 0
    % errRot = real(decision)*Q - imag(decision)*I;                                   % Phase detector

    errQuad = sign(Q)*(I - real(decision));    % This is sign(Q)(I-^I)

    acqDetectBuffer = [acqDetectBuffer; abs(errRot) + abs(errQuad)];

    theta = theta + k * errRot;
    beta = beta + k * errQuad;

    % Save arrays for plotting results
    errsRot = [errsRot; errRot];
    errsQuad = [errsQuad; errQuad];
    thetas = [thetas; theta];
    betas = [betas; beta];
    Iout = [Iout; I];
    Qout = [Qout; Q];
end

%% Plot the starting constellation
figure;
plot(real(sampsNorm + noise), imag(sampsNorm + noise), 'r.');
axis equal;
grid on;
title('Starting Constellation');

%% Plot the convergence of phase rotation in the loop
figure;
plot(thetas);
title('Convergence of Phase Rotation in Loop');
xlabel('Sample Number');
grid on;

%% Plot the convergence of beta (quadrature error) in the loop
figure;
plot(betas);
title('Convergence of Beta (Quadrature Error) in Loop');
xlabel('Sample Number');
grid on;

%% Plot the corrected constellation
figure;
plot(Iout, Qout, 'r.');
axis equal;
grid on;
title('Corrected Constellation');

%% Plot the error signal for acquisition detection

figure;

% This is the combined absolute error from rotation and quadrature error which can
% be used to determine "acquired" once error is below a set threshold:
acqDetectBuff = abs(errsRot) + abs(errsQuad);

plot(acqDetectBuff, 'DisplayName', 'Sum of |Rot| and |Quad| Errors');
hold on
plot(filter(ones(1, acqsamps), acqsamps, acqDetectBuff), 'DisplayName', 'Moving Avg');
title('Error Signal for Acquisition Detection');
xlabel('Sample Number');
ylabel('Error');
legend();





%% Quadrature Correction Demo on QAM Waveform
% Create random QAM data symbols with a static phase offset to demonstrate

num_symb = 2000;
symbols = qam(num_symb, 16); % Assuming you have a function 'qam' that generates QAM symbols

% Add a phase offset
phase = -0.8;

% Add quadrature error
quad = -0.2;

I = real(symbols);
Q = imag(symbols);

tx_nom = I + 1i * Q;

tx = ((I + Q * sin(quad)) + 1i * Q) * exp(1i * phase);
fprintf('STD of QAM waveform: %.2f\n', std(tx));

figure;
plot(real(tx_nom), imag(tx_nom), 'r+');
title('16QAM Constellation');
axis equal;

figure;
plot(real(tx), imag(tx), 'o');
title('16QAM with Static Phase Offset and Quadrature Error');
axis equal;

%% Normalize waveform
txNorm = std(symbols) / std(tx) * tx;

decisionI = (round(real(txNorm) / 2 + 0.5) - 0.5) * 2;
decisionQ = (round(imag(txNorm) / 2 + 0.5) - 0.5) * 2;
decisions = decisionI + 1i * decisionQ;

figure;
plot(real(txNorm), imag(txNorm), 'o', 'DisplayName', 'Rx Symbols');
hold on;
plot(decisionI, decisionQ, '+', 'DisplayName', 'Closest Decisions');
axis equal;
title('16QAM with Static Phase Offset and Quadrature Error');
legend();

%% Alternate phase estimator using quadrature samples and sign of real
rotCorrect = sign(real(tx)) .* (imag(tx) - imag(decisions));

%% Plot the result of the alternate phase estimator
figure;
plot(rotCorrect);
fprintf('Mean of rotCorrect: %.4f\n', mean(rotCorrect));
fprintf('True phase value: %.4f\n', phase);

%% Correction loop to remove rotation (Type 1 first-order loop)

accumQ = 0;
k = 0.5;

% Save results of loop simulation for plotting
errs = [];
derotatedSamples = [];
accums = [];

for count = 1:length(txNorm)
    if count > 30 % Gear shift after acquisition
        k = 0.001;
    end
    
    sample_rotated = txNorm(count) * exp(-1j * accumQ);
    derotatedSamples = [derotatedSamples; sample_rotated];
    
    decisionI = round(real(sample_rotated) / 2 + 0.5) * 2 - 0.5;
    decisionQ = round(imag(sample_rotated) / 2 + 0.5) * 2 - 0.5;
    
    % New phase detector that will drive Q difference to 0 as desired for subsequent quadrature removal
    % Can remove scaling and phase noise increases
    errQ = sign(real(sample_rotated)) * (imag(sample_rotated) - decisionQ) / abs(decisionI + 1j * decisionQ);
    
    % Phase detector given by ^I*Q - I*^Q will leave a rotation offset due to quadrature imbalance
    % Although this may still be better for higher order QAM due to consistent error result for decisions closer to 
    % imaginary axis, and when combined with quadrature correction in loop should still converge. 
    % errQ = real(decision) * imag(sample_rotated) - imag(decision) * real(sample_rotated);
    
    errs = [errs; errQ];
    accumQ = accumQ + k * errQ;
    accums = [accums; accumQ];
end

%% Plot the result of the correction loop
figure;
plot(real(derotatedSamples(100:end)), imag(derotatedSamples(100:end)), 'bo', 'DisplayName', 'Rx Symbols');
hold on;
plot(real(tx_nom), imag(tx_nom), 'r+', 'DisplayName', 'Closest Decisions');
axis equal;
title('16QAM with Phase Offset Corrected (IQ Error Remains)');
legend();

%% Plot the measured quadrature error from the correction loop
figure;
plot(accums);
title('Measured Quadrature Error from Correction Loop');
xlabel('Sample');
ylabel('Quadrature Error');
