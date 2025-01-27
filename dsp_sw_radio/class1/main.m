%% expample using spPoly
spPoly([1, pi], [1,2,1], 'z', 2);

%% Accumulator Step Response

% Create time vector
t = 0:0.1:2.4;

accum = tf(0.1*[1, 0], [1, -1],.1);

% Compute the step response
y = step(accum, t);

% Plot the step response
figure;
plot(t, y, 'ro');
grid on;
title('Step Response of Accumulator with fs = 10 Hz');
xlabel('Time (seconds)');
ylabel('Amplitude');

%% Magnitude and Phase Response
% Frequency Response of 1 + z^(-1)

% Compute the frequency response
 freqz([1, 1], 1, 'whole')

% Plot the magnitude response
figure;
subplot(2,1,1);
plot(w, abs(H)); 
ylabel('Magnitude');
title('Frequency Response of 1 + z^{-1}');
grid on;

% Plot the phase response
subplot(2,1,2);
plot(w, angle(H) * 180/pi);
xlabel('Normalized Angular Frequency (2\pi = Fs)');
ylabel('Phase [deg]');
grid on;

%% Laplace Transform 
% This uses the example on slide 35 in the class 1 presentation to further 
% show the impulse response, Laplace Transform and Frequency response for the filter given.
syms s
t = linspace(0, 60, 1000); % adjust the time vector as needed
impulse_response = ilaplace(1/((s+0.2+0.5i)*(s+0.2-0.5i)));
response_values = subs(impulse_response, t);
figure;
plot(t, response_values);
xlabel('sec');
title('Impulse Response for x(t)');
grid on;

%% 
num = [1];
den = conv([1, 0.2+0.5i], [1, 0.2-0.5i]);

sys = tf(num, den);

t = linspace(0, 60, 512);
y = impulse(sys, t);

figure;
plot(t, y);
title('Impulse Response - Continuous System');
xlabel('Time [s]');
ylabel('Magnitude');
grid on;

%%
num = [1];
den = conv([1, 0.2+0.5i], [1, 0.2-0.5i]);

[h,w] = freqs(num, den, linspace(-10, 10, 500));
figure;
plot(w, abs(h));
xlabel('Frequency [rad/sec]');
ylabel('Magnitude');
title('Frequency Response (Mag) for H(s)');
axis([-1, 1, 0, 6]);    % scaled to match what is shown in class presentation slide 35
grid on;

%% QAM-16 Constellation
qam = qammod(randi([0 15],100,1),16,"gray","PlotConstellation",true);


