clear all
close all
clc;

M=16; % 16-QAM

% Create random QAM data symbols with a static phase offset
numSymbols = 200;

% Create a random data source
data = randi([0 15], 1, numSymbols);
phase = pi/20;
% Map symbols using the constellation
symbols = qammod(data, M);

% Display the constellation diagram
scatterplot(symbols);
title('16-QAM with Static Phase Offset');
axis equal;

% Add a static phase offset
phaseOffset = pi/20;
disp(['Phase error introduced is ' num2str(phaseOffset) ' radians']);
tx = symbols .* exp(1i * phaseOffset);

%% Decision Directed Phase Detector

% Generator decisions:
% Round I and Q to be closest to +/-1 and +/-3
% This is done by first dividing by 2 to have
% closest to +/-1.5 and +/-0.5
% Then shifting by 0.5 to get 1, 0, 1, 2 in order
% to use the integer round function, and then once
% rounded, return to the original levels:

% Create a complex signal tx (16-QAM with static phase offset)
tx = symbols .* exp(1i * phaseOffset);

% Generator decisions
decisionI = (round(real(tx)/2 + 0.5) - 0.5) * 2;
decisionQ = (round(imag(tx)/2 + 0.5) - 0.5) * 2;
decisions = decisionI + 1i * decisionQ;

% Plot the decisions
figure;
plot(decisionI, decisionQ, '.');
hold on;
plot(real(decisions(1)), imag(decisions(1)), '+');
hold on
plot(real(tx(1)), imag(tx(1)), '+');
title('16-QAM Decisions');
axis equal;

% Calculate average magnitude of decisions and tx
A = mean(abs(decisions));
disp(['Average magnitude of decisions is ' num2str(A)]);
disp(['Average magnitude of tx is ' num2str(mean(abs(tx)))]);



%% Compute the phase detector error between actual samples and closest decision:

% Showing the two approaches of the same computation, the first
% demonstrates processing on the I and Q samples as would be done in 
% implementation, the second shows the direct calculation on the 
% complex samples, which as demonstrated by tic/toc is faster for 
% simulation; both provide the same result.

% Approach 1:
tic;
ph_error = decisionI .* imag(conj(tx) .* tx) - decisionQ .* real(conj(tx) .* tx);
toc;

% Approach 2:
tic;
ph_error = imag(conj(decisions) .* tx);
toc;

disp(['True phase error = ' num2str(phase)]);
disp(['Proportional phase error measured = ' num2str(mean(ph_error))]);
disp(['Detector gain = ' num2str(mean(ph_error) / phase)]);
disp(['Proportional phase error predicted by A^2 sin(theta) = ' num2str(mean(abs(decisions).^2) * sin(phase))]);

figure;
plot(ph_error);
title('Sample by Sample Computed Phase Error');
xlabel('Sample Number');
ylabel('Proportional Phase Error');


%% FFT of phase error result

npts = 2^14;

% Note: cannot use small angle approximation and just take FFT of ph_error directly!
ph_noise = fftshift(fft(exp(1i * ph_error) / length(ph_error), npts));
frange = (-npts/2):(npts/2-1);

figure;
plot(2 * pi * frange / npts, 20 * log10(abs(ph_noise)));
grid on;
ylabel('dB');
xlabel('Normalized Frequency (rad/sample)');
title('Self-Noise Spectrum');


%% Magnitude-scaled phase measurement

disp(['True phase error = ', num2str(phase, '%0.4f')]);

ph_error = imag(conj(decisions) .* tx) ./ (abs(decisions).^2)

disp(['Estimated phase error = ', num2str(mean(ph_error), '%0.4f')]);

figure;
plot(ph_error);
title('Sample by Sample Computed Phase Error with Mag Scaling');
xlabel('Sample Number');
ylabel('Phase Error');

%% Same with a frequency offset

% Create random data source
num_symb = 200;

data = randi([0, 15], 1, num_symb);
symbols = qammod(data, 16);

% Increase phase increment on every count
phase_step = pi/200;
phase_ramp = phase_step*(1:num_symb);
tx = symbols .* exp(1j*phase_ramp);

% Round I and Q to be closest to +/-1 and +/-3
decisionI = (round(real(symbols)/2 + 0.5) - 0.5) * 2;
decisionQ = (round(imag(symbols)/2 + 0.5) - 0.5) * 2;
decisions = decisionI + 1i*decisionQ;

% Phasor error computation
ph_error = imag(conj(decisions) .* tx) ./ (abs(decisions).^2);

figure;
plot(ph_error, 'DisplayName', 'Estimated Error');
hold on;
plot(phase_ramp, 'DisplayName', 'True Error');
hold off;
legend;
title('Sample by Sample Computed Phase Error');
xlabel('Sample Number');
ylabel('Phase Error');

%% Demonstrating carrier offset measurement multiple samples/symbol with pulse shaping

% First upsample the data
% with a raised cosine pulse shaping filter (see Class 1):
num_symb = 24000;
sym_rate = 1e3;
sym_length = 26;      % Duration of impulse response
oversamp = 4;         % Samples per symbol
fs = sym_rate * oversamp;   % Sampling rate
alpha = 0.5;          % Roll-off factor (set between 0 and 1)
data = randi([0, 15], 1, num_symb);
symbols = qammod(data, 16);

h_rrc = rcosdesign(alpha, sym_length, oversamp, 'sqrt')/sqrt(oversamp);
h_rc =conv(h_rrc,h_rrc);
[h, w] =freqz(h_rc, 1, "half",fs);
index = linspace(-sym_length/2, sym_length/2, sym_length * oversamp);

tx = zeros(1, num_symb * oversamp);
tx(1:oversamp:end) = symbols;

% Pass data through tx pulse shaping filter:
tx_shaped = filter(h_rc, 1, tx);

%% IQ Diagram of Modulated Waveform with No Carrier Offset

figure;
plot(real(tx_shaped), imag(tx_shaped));
axis equal;
title('IQ Diagram of Modulated Waveform with No Carrier Offset');

%% Discriminator Output, No Frequency Offset

product = tx_shaped .* conj(circshift(tx_shaped, [0, 1]));

fout = imag(product);
disp(['The average out of the discriminator is ', num2str(mean(fout), '%0.3f')]);
disp(['The std out of the discriminator is ', num2str(std(fout), '%0.3f')]);

figure;
plot(fout);
hold on;
plot([0, length(tx_shaped)], [mean(fout), mean(fout)]);
xlabel('Time (samples)');
title('Case 1: Discriminator Output, No Frequency Offset');

%% Induce a Frequency Offset

foffset = sym_rate * 0.1;
t = (0:length(tx_shaped)-1)/fs;

tx_offset = tx_shaped .* exp(1j*2*pi*foffset*t);

figure;
plot(real(tx_offset), imag(tx_offset));
axis equal;
title('IQ Diagram of Modulated Waveform with Carrier Offset');

%% Calculate Discriminator Output with Frequency Offset

fout = imag(tx_offset .* conj(circshift(tx_offset, [0, 1])));
disp(['The average out of the discriminator is ', num2str(mean(fout), '%0.3f')]);
disp(['The std out of the discriminator is ', num2str(std(fout), '%0.3f')]);

figure;
plot(fout);
hold on;
plot([0, length(tx_offset)], [mean(fout), mean(fout)]);
title('Case 1: Discriminator Output, Frequency Offset');

%% FFT of Discriminator Output (Frequency Noise)

npts = 2^14;
freq_noise = fftshift(fft(fout / length(fout), npts));
frange = -npts/2 : npts/2 - 1;

figure;
plot(2*pi*frange/npts, 20*log10(abs(freq_noise)));
grid on;
ylabel('dB');
xlabel('Normalized Frequency (rad/sample)');
title('Discriminator Spectrum');


%% Timing Error Detection: Eye Diagram

num_cycles = 4;      % number of symbols to display in eye diagram
windows = 200;       % a window is one path across display shown 
upsample = 16;

% Resample data to 64x per symbol to emulate continuous waveform
tx_resamp = resample(tx_shaped, upsample, 1);

samp_per_win = oversamp * upsample * num_cycles;

% Divide by number of samples per win and then 
% pad zeros to next higher multiple using tx_eye = reshape(tx_shaped, N, [])
N = floor(length(tx_resamp) / samp_per_win);
tx_eye = [tx_resamp(:); zeros(N * samp_per_win - length(tx_resamp), 1)];
grouped = reshape(tx_eye, samp_per_win, N);

transient = sym_length / 2;
eye = real(grouped);

% Create an x-axis in samples np.shape(eye) gives the
% 2-dimensional size of the eye data and the first element
% is the interpolated number of samples along the x-axis
nsamps = size(eye, 1);
xaxis = (1:nsamps) / upsample;

figure;

% Plot showing continuous trajectory of 
plot(xaxis, eye(:, transient : transient + windows));

% Actual sample locations
hold on;
plot(xaxis(1 : upsample : end), eye(1 : upsample : end, transient : transient + windows), 'b.');
hold off;

title('Eye Diagram');
xlabel('Samples');
grid on;

%% Algorithms for timing recovery
% Gardner 

% Generate all offsets over symbol duration
offsets = 0:upsample*oversamp-1;

% Rotate offsets to center error detection plot
offsets = circshift(offsets, -25);

% Increment through each offset and compute the average timing error
ted_results = gardner_ted(tx_resamp, upsample*oversamp, offsets);

figure;
plot(-31:32, ted_results);
title('Gardner TED Results for QAM Waveform');
xlabel('Sample Offset');
ylabel('Measured Timing Error');
grid on;


% Mueller and Mueller 

% Increment through each offset and compute the average timing error
mm_results = mandm_ted(tx_resamp, upsample*oversamp, offsets);

figure;
plot(-31:32, mm_results);
title('M&M Results for QAM Waveform');
xlabel('Sample Offset');
ylabel('Measured Timing Error');
grid on;

%% Carrier offset measurement under timing offset conditions
% consider different conditions of time offset by changing n1 below and running the subsequent cells
n1 = 0; % change n1=0,1,2,3 to select other offsets with reference to the eye diagram above
tx_time_offset = 4 * tx_shaped(n1+1 : 4 : end); % scaled to normalize with decision boundaries

figure;
plot(real(tx_time_offset), imag(tx_time_offset), 'o');
title('16QAM w/ Timing Offset');
axis equal;

% add phase offset:
phase_offset = 0.51;
tx_ph_time_offset = tx_time_offset * exp(1j * phase_offset);
tx_ph = tx_shaped * exp(1j * phase_offset);

figure;
plot(real(tx_ph_time_offset), imag(tx_ph_time_offset), 'o');
title('16QAM w/ Timing and Phase Offset');
axis equal;

% Decision Directed Phase Detector result (using magnitude scaled implementation introduced above)
decisionI = (round(real(tx_ph_time_offset)/2 + 0.5) - 0.5) * 2;
decisionQ = (round(imag(tx_ph_time_offset)/2 + 0.5) - 0.5) * 2;
decisions = decisionI + 1i * decisionQ;

ph_error = imag(conj(decisions) .* tx_ph_time_offset) ./ (abs(decisions).^2);

fprintf('The average out of the discriminator is %0.3f\n', mean(ph_error));
fprintf('The std out of the discriminator is %0.3f\n', std(ph_error));

figure;
plot(ph_error);
hold on;
plot([0, length(ph_error)], [mean(ph_error), mean(ph_error)]);
title('Sample by Sample and Average Out of Phase Detector');

figure;
plot((-length(ph_error)/2:length(ph_error)/2-1), fftshift(20*log10(abs(fft(ph_error)))));
title('FFT of Phase Detector Output');
xlabel('Frequency Bin');
ylabel('dB');
axis([-length(ph_error)/2, length(ph_error)/2, 10, 70]);


%% demonstrating x^4 recovery for QAM waveform above

% the phase of the resulting detector output will be 4x the actual phase error, and
% similarly the frequency of the resulting detector output will be 4x the actual frequency error.
% For the tracking loop, we really just need an error metric, so the factor of 4 becomes part of
% the loop gain.

% Below as the FFT shows the strong DC component, as expected for the case of a static phase error.
% If this phase was changing with time (frequency), we would see the resulting shift in frequency 
% in the FFT.

error = (tx_ph .* exp(1j*2*pi*0.01*(1:length(tx_ph)))).^4;

figure;

subplot(2,1,1);
plot((-length(tx_ph)/2:length(tx_ph)/2-1), fftshift(20*log10(abs(fft(tx_ph .* kaiser(length(tx_ph),8).')))));
grid on;
axis([-length(error)/2, length(error)/2, -80, 60]);
title('FFT of 16-QAM Waveform');
ylabel('dB');
xlabel('Frequency Bin');

subplot(2,1,2);
plot((-length(error)/2:length(error)/2-1), fftshift(20*log10(abs(fft(error)))));
title('FFT of (QAM)^4 Output');
ylabel('dB');
xlabel('Frequency Bin');
axis([-length(error)/2, length(error)/2, 20, 90]);
grid on;

sgtitle('Demonstrating x^4 Recovery for QAM Waveform');

%% complex IQ plot of error signal
fprintf('%.4f %.4f\n', mean(real(error)), mean(imag(error)));
figure;
plot(real(error), imag(error), 'LineWidth', 0.5);
hold on
arrowScale = 4 * mean(real(error)) / abs(mean(error));
arrowDirection = 4 * mean(imag(error)) / abs(mean(error));
quiver(0, 0, arrowScale, arrowDirection, 'Color', 'red', 'LineWidth', 2);
title('(16QAM)^4 w/ Timing and Phase Offset');
axis equal;
grid on;





%% Carrier Tracking Loop
% 𝐺ol(𝑧)=2𝜋𝐴^2(𝑇+τ2)/τ1 *(𝑧−(τ2/(Τ+𝜏2))/(𝑧−1)^2

% Loop Parameters
A = 1;                  % magnitude of complex input
T = 1/12e6;             % symbol duration in seconds (inverse of loop update rate)
w3db = 2*pi*76.32e3;    % 3 dB loop BW
zeta = 0.5;             % damping factor

Kv = 2*pi/T;            % NCO gain
Kd = A^2;               % Discriminator gain
alpha = 1 - 2*zeta^2;
wn = w3db / sqrt(alpha + sqrt(alpha^2 + 1));  % natural frequency
tau1 = Kv * Kd / wn^2;  % loop filter time constant
tau2 = 2 * zeta / wn;    % Loop filter time constant

% The numerator z in NCO cancels with denominator z from parasitic delay in loop
H_nco = tf([Kv*T], [1, -1], T);  % NCO transfer function (phase accumulator) without numerator z
H_lf = tf([(T + tau2) / tau1, -tau2 / tau1], [1, -1], T);  % PI Loop filter

Gol = H_nco * H_lf * Kd;
disp(Gol);

% Plot Bode Plot
figure;
bode(Gol);

% Print open-loop zero and open-loop gain
fprintf('Open loop zero = %.4f\n', tau2 / (T + tau2));
fprintf('Open loop gain = %.4f\n', 2 * pi * A^2 * (T + tau2) / tau1);

% Print additional loop parameters
fprintf('Natural frequency = %.4f rad/sec\n', wn / (2 * pi));
fprintf('Loop bandwidth = %.4f rad/sec\n', w3db);
fprintf('T/tau1 = %.2f us\n', T / tau1 * 1e6);
fprintf('tau1/tau2 = %.4f\n', tau1 / tau2);
fprintf('tau1 = %.2f us\n', tau1 * 1e6);
fprintf('tau2 = %.2f us\n', tau2 * 1e6);
disp(H_nco);
disp(Kv * Kd / (2 * pi));
disp(T - tau2);
disp(tau2 / tau1);
disp(T / tau1);
disp(tau2 / T);

% Step response
figure;
t = (0:349) * T;
yout = step(1 / (1 + Gol), t);
plot(t / T, yout);
grid on;
xlabel('Symbol Number');
ylabel('Phase Error [radians]');
title('Step Response');

% Ramp response
integrator = tf([1, 0], [1, -1], T);  % transfer function of z/(z-1)
mag = 2 * pi * 55e3 / 12e6;  % magnitude of step based on freq offset, normalized to loop update rate

% Loop gain 0.0324 as given above:
[yout, t] = step(mag * integrator / (1 + Gol), (0:349) * T);

% Compare to loop gain 0.076 as shown in the presentation
% I do this by dividing Gol by 0.0324 and then multiplying it by 0.076 in the closed
% loop transfer function input below:
[yout2, t2] = step(mag * integrator / (1 + (0.076 / 0.0324) * Gol), (0:349) * T);

figure;
plot(t / T, yout, 'DisplayName', 'Gain = 0.0324');
hold on;
plot(t2 / T, yout2, 'DisplayName', 'Gain = 0.0760');
grid on;
xlabel('Symbol Number');
ylabel('Phase Error [radians]');
title('Ramp Response');
legend;
axis([0, 350, -0.2, 0.8]);

% Pole-zero map
sys_cl = 1 / (1 + Gol);
figure;
pzmap(sys_cl);
axis([0.95, 1.02, -0.05, 0.05]);
grid on;


function ted_results = gardner_ted(tx, n, offsets)
    % tx: oversampled complex waveform
    % n: oversampling rate
    % offsets: array of sample offsets to test
    
    ted_results = zeros(size(offsets));
    
    for i = 1:length(offsets)
        offset = offsets(i);
        
        % Downsample to 2 samples per symbol with timing offset
        tx2 = tx(offset + 1 : n/2 : end);
        
        % Generate a prompt, late, and early each offset by 1 sample
        late = tx2(3:end);
        early = tx2(1:end-2);
        prompt = tx2(2:end-1);
        
        % Compute and store the Gardner Error result
        tmp = (late - early);
        ted_results(i) = mean(real(conj(prompt(1:2:end)).*tmp(1:2:end)));
    end
end

function mm_results = mandm_ted(tx, n, offsets)
    % tx: oversampled waveform
    % n: oversampling rate
    % offsets: array of sample offsets to test
    
    mm_results = zeros(size(offsets));
    
    for i = 1:length(offsets)
        offset = offsets(i);
        
        % Downsample to 1 sample per symbol with timing offset
        tx2 = tx(offset + 1 : n : end);
        
        % Compute M&M on real axis
        sign_tx2_real = sign(real(tx2));
        mm_real = real(tx2(2:end)).*sign_tx2_real(1:end-1) - real(tx2(1:end-1)).*sign_tx2_real(2:end);
        
        % Compute M&M on imag axis
        sign_tx2_imag = sign(imag(tx2));
        mm_imag = imag(tx2(2:end)).*sign_tx2_imag(1:end-1) - imag(tx2(1:end-1)).*sign_tx2_imag(2:end);
        
        % Store the M&M results
        mm_results(i) = mean(mm_real + mm_imag);
    end
end