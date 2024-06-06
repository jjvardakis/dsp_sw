% Loop Parameters
A = 1;                % Magnitude of complex input
T = 1;                % Symbol duration in seconds (inverse of loop update rate)
I = 0.072;            % Fast: 0.075 (5 symbols)   Slow: 
P = 0.21;             % Fast: 0.21 (5 symbols) Slow

% Calculate the open-loop transfer function
Gol_wide = ctrack(A, T, I, P);
disp(Gol_wide);

% Display open loop zero and gain
disp(['Open loop zero = ', num2str((P - I * T) / P, '%0.4f')]);
disp(['Open loop gain = ', num2str(2 * pi * P * A^2, '%0.4f')]);


%% Define the system transfer function
Gol_wide = ctrack(A, T, I, P);

% Plot the step response
figure;
[yout,t] = step(1 / (1 + Gol_wide), (0:19) * T);
plot(t / T, yout);
grid on;
xlabel('Symbol Number');
ylabel('Phase Error [radians]');
title('Step Response - Wide Loop BW');
axis([0, 20, -0.5, 1.3]);

%% Define the system transfer function
Gol_wide = ctrack(A, T, I, P);

% Compute the closed-loop transfer function
Gcl_wide = minreal(1 / (1 + Gol_wide));

% Plot the poles and zeros
figure;
pzmap(Gcl_wide);
grid on;

%% Define the frequency range for the Bode plot
omega = logspace(-2, log10(pi), 512);

% Compute the magnitude and phase of the closed-loop transfer function
% [mag, phase] = bode(Gcl_wide, omega);
bode(Gcl_wide, omega)

% % Plot the Bode plot
% figure;
% semilogx(omega, 20*log10(mag));
grid on;
xlabel('Normalized Radian Frequency (cycles/sample)');
ylabel('Magnitude (dB)');
title('Closed Loop Frequency Response - Wide Loop BW');

%% Compute the frequency response using freqz
[h, w] = freqz(numden(Gcl_wide), 1);

% Plot the frequency response
figure;
semilogx(w, 20*log10(abs(h)));
grid on;
xlabel('Normalized Radian Frequency (cycles/sample)');
ylabel('Magnitude (dB)');
title('Closed Loop Frequency Response - Wide Loop BW');
axis(axis);

%% Root Locus: showing us the possible pole locations as we vary gain:
% "gain" is a constant K in the closed loop system such that the closed loop
% transfer function is 1/(1+K*Gol)
figure;
rlocus(Gol_wide);
axis([0 1 -0.5 0.5]);



%% Medium Loop BW
% Medium loop
% Loop Parameters
A = 1;                % magnitude of complex input
T = 1;                % symbol duration in seconds (inverse of loop update rate)
I = 0.010;            % Fast: 0.072 (5 symbols)   Slow: 0.0001  (200 symbols)
P = 0.07;             % Fast: 0.21    (5 symbols)    Slow: 0.007   (200 symbols)

Gol_mid = ctrack(A, T, I, P);
fprintf('Open loop zero = %.4f\n', (P-I*T)/P);
fprintf('Open loop gain = %.4f\n', 2*pi*P*A^2);

figure;
[yout,t] = step(1/(1+Gol_mid), (0:49)*T);
plot(t/T, yout);
grid on;
xlabel('Symbol Number');
ylabel('Phase Error [radians]');
title('Step Response - Medium Loop BW');
axis([0, 50, -0.5, 1.3]);

%% Calculate closed-loop transfer function
Gcl_mid = minreal(1/(1+Gol_mid));
disp(Gcl_mid);

% Plot poles and zeros of closed-loop system
figure;
pzmap(Gcl_mid);
grid on

%% Plot closed-loop frequency response
figure;
bode(Gcl_mid, {0.01, pi});
grid on;
title('Closed Loop Frequency Repsonse - Medium Loop BW');

%% Narrow Loop BW
% Define the loop parameters
A = 1;            % magnitude of complex input
T = 1;            % symbol duration in seconds (inverse of loop update rate)
I = 0.0001;       % Fast: 0.072 (5 symbols)   Slow: 0.0001  (200 symbols)
P = 0.007;        % Fast: .21    (5 symbols)   Slow: 0.007   (200 symbols)

% Create the transfer function for the open-loop system
Gol_narrow = ctrack(A, T, I, P);

% Display open-loop zero and gain
fprintf('Open loop zero = %.4f\n', (P - I * T) / P);
fprintf('Open loop gain = %.4f\n', 2 * pi * P * A^2);

% Plot the step response of the closed-loop system
figure;
t = 0:T:349 * T;
yout = step(1 / (1 + Gol_narrow), t);
plot(t / T, yout);
grid on;
xlabel('Symbol Number');
ylabel('Phase Error [radians]');
title('Step Response - Narrow Loop BW');
axis([0, 350, -0.5, 1.3]);

% Obtain the closed-loop tr
% Poles and zeros of closed-loop system
Gcl_narrow = minreal(1 / (1 + 2 * Gol_narrow));

% Display poles and zeros
figure;
pzmap(Gcl_narrow);
grid on;
title('Poles and Zeros - Narrow Loop BW');

%% Calculate the frequency response
bode(Gcl_narrow);
grid on;
xlabel('Normalized Radian Frequency (cycles/sample)');
ylabel('Magnitude (dB)');
title('Closed Loop Frequency Response - Narrow Loop BW');


%% Cancelling pole/zero pairs:minreal
Gcl = Gol_narrow/(1+Gol_narrow)
disp(Gcl)

% Calculate the closed-loop transfer function and cancel common factors
Gcl = minreal(Gol_narrow/(1 + Gol_narrow));

% Display the result
disp(Gcl);

% Plot the poles and zeros
figure;
pzmap(Gcl);
grid on;

%% Plot the root locus
figure;
rlocus(Gol_narrow);
grid on;

%% Plot the step response
figure;
t = 0:T:349*T;
yout = step(Gcl, t);
plot(t/T, yout);
grid on;
xlabel('Symbol Number');
ylabel('Phase Error [radians]');
title('Step Response - Narrow Loop BW');
axis([0 350 -0.5 1.3]);

%% Plot the closed-loop frequency response
figure;
bode(Gcl_narrow, logspace(-2, log10(pi), 512));
grid on;
xlabel('Normalized Radian Frequency (cycles/sample)');
ylabel('Magnitude (dB)');
title('Closed Loop Frequency Response - Narrow Loop BW');


%% Noise optimization
% Create phase noise samples

N = 2^15;
Fs = 1.0;
Sin = ones(N, 1) + 1j * 0.001 * ones(N, 1);

% Define frequency points and corresponding power levels
fpoint = [1e-2, 1e-1, 0.49];
npower = [-15, -45, -59];

% Generate phase noise
[F, P, pnoise] = ph_noise(Sin, Fs, fpoint, npower);

figure()
plot(pnoise)
ylabel('Phase error (radians)')
xlabel('Sample Number')
title('Added Phase Noise')

% Plot the phase noise
figure;
plot(F, P);
grid on;
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Phase Noise');

%% Apply the phase noise correction
% Add AWGN noise
% npower = -40; % Noise power in dBc
% ndensity = npower - 10*log10(Fs);
% 
% awgn = 10^((npower-3)/20) * (randn(N, 1) + 1i * randn(N, 1));
% 
% noise = pnoise + awgn;
% 
% figure;
% 
% % Calculate input length
% N = numel(noise);
% 
% % M is the frequency resolution
% if mod(N, 2) % N Odd
%     M = (N + 1) / 2 + 1;
% else         % N Even
%     M = N / 2 + 1;
% end
% F = linspace(0, Fs/2, M);
% dF = (F(end) - F(end-1)) * ones(1, M); % Delta F (Multiplied by M to create 1xM array)
% X1 = fft(noise);
% 
% semilogx(F(2:end), 10*log10((abs(X1(1:M-1))/max(abs(X1(1:M-1)))).^2 / dF(1)));
% grid on;
% title('Spectral Density of Carrier Noise (Phase Noise + AWGN)');
% xlabel('Frequency Offset (cycles/sample)');
% ylabel('Rel. Power (dBc)');
% % axis([1e-2, 0.5, -70, 0]);
