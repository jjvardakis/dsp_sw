%% Step Response - 2nd Order System
zetas = [0.2, 0.7, 1.5];
t = linspace(0, 30, 1000);

figure;
for zeta = zetas
    sys = tf(1, [1, 2 * zeta, 1]);
    [yout,t ] = step(sys, t);
    plot(t, yout, 'DisplayName', sprintf('zeta = %.1f', zeta));
    hold on;
end

legend;
xlabel('Time [s]');
ylabel('Amplitude');
title('Step Response 2nd Order System');

%% Roots of the characteristic equation
roots([1, 2, 2, 0]) % finding the poles of the open loop gain (the roots of the denominator)


%% Nyquist Plot Example
K = 5;
gol = tf(K, [1, 2, 2, 0]);

figure;
nyquist(gol);
axis([-K/2 - 1, K/2 + 1, -K/2 - 1, K/2 + 1]);

%% Step response
K = 3.6;
gol = tf(K, [1, 2, 2, 0]);

figure;
T = linspace(0, 120, 500);
[yout, n] = step(gol/(1 + gol), T);
plot(n, yout, 'LineWidth', 1.5, 'DisplayName', ['K = ', num2str(K)]);
title('Step Response');
xlabel('Time (sec)');
legend();
grid on;

%% Using a Step Response test to estimate loop BW
wn = 2 * pi * 1000;
zeta1 = 0.5;
zeta2 = 0.7;
zeta3 = 1.0;

gol1 = tf(wn^2, [1, 2 * zeta1 * wn, wn^2]);
gol2 = tf(wn^2, [1, 2 * zeta2 * wn, wn^2]);
gol3 = tf(wn^2, [1, 2 * zeta3 * wn, wn^2]);

T = linspace(0, 0.002, 500);
[yout1, t1] = step(2 * gol1/(1 + gol1), T);   
[yout2, t2] = step(2 * gol2/(1 + gol2), T);
[yout3, t3] = step(2 * gol3/(1 + gol3), T);

figure;
plot(t1 * 1000, yout1, 'LineWidth', 1.5, 'DisplayName', '\zeta = 0.5');
hold on;
plot(t2 * 1000, yout2, 'LineWidth', 1.5, 'DisplayName', '\zeta = 0.7');
plot(t3 * 1000, yout3, 'LineWidth', 1.5, 'DisplayName', '\zeta = 1.0');
plot(t1 * 1000, 0.9 * ones(size(t1)), 'r--');
plot(t1 * 1000, 0.1 * ones(size(t1)), 'r--');
title('Step Response of 2nd order System');
xlabel('Time (mS)');
ylabel('Normalized Magnitude');
legend();
grid on;

figure;
[mag1, ph1, freq1] = bode(gol1/(1 + gol1));   
[mag2, ph2, freq2] = bode(gol2/(1 + gol2));
[mag3, ph3, freq3] = bode(gol3/(1 + gol3));
mag1 = squeeze(mag1);
mag2 = squeeze(mag2);
mag3 = squeeze(mag3);
subplot(2,1,1);
semilogx(freq1/(2*pi), db(mag1) - db(mag1(1)), 'LineWidth', 1.5, 'DisplayName', '\zeta = 0.5');
hold on;
semilogx(freq2/(2*pi), db(mag2) - db(mag2(1)), 'LineWidth', 1.5, 'DisplayName', '\zeta = 0.7');
semilogx(freq3/(2*pi), db(mag3) - db(mag3(1)), 'LineWidth', 1.5, 'DisplayName', '\zeta = 1.0');
title('Frequency Response of 2nd order System');
axis([ 0 16000 -48 5])
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend();
grid on;

subplot(2,1,2);
plot(freq1/(2*pi), db(mag1) - db(mag1(1)), 'LineWidth', 1.5, 'DisplayName', '\zeta = 0.5');
hold on;
plot(freq2/(2*pi), db(mag2) - db(mag2(1)), 'LineWidth', 1.5, 'DisplayName', '\zeta = 0.7');
plot(freq3/(2*pi), db(mag3) - db(mag3(1)), 'LineWidth', 1.5, 'DisplayName', '\zeta = 1.0');
plot(freq3/(2*pi), -3 * ones(size(freq3)), 'r--');
axis([ 0 16000 -48 5])
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend();
grid on;

%% Ramp response
K = 1;
gol = tf(K, [1, 2, 2, 0]);

integrator = tf(1, [1, 0]);

figure;
T = linspace(0, 20, 500);
[yout, n] = step(integrator * gol/(1 + gol), T);
plot(n, yout, 'LineWidth', 1.5, 'DisplayName', ['K = ' num2str(K)]);
title('Step Response');
xlabel('Time (sec)');
legend();
grid on;

%% Root locus
K = 1;
gol = tf(K, [1, 2, 2, 0]);

figure;
rlocus(gol);
title('Root Locus');

%% Pole Zero Map
K = 1;
gol = tf(K, [1, 2, 2, 0]);

figure;
pzmap(1/(1+K*gol),'r')
grid on

%% Mapping
T = pi/10;
gol_d = c2d(gol, T, 'impulse');
disp(gol_d);
figure;
pzmap(gol_d);

% Comparing frequency response of the continuous and discrete systems
figure;
bode(gol_d)
hold on
bode(gol)

p2 = exp((-1-1j)*T);
p3 = exp((-1+1j)*T);
disp(p2);
disp(p3);
conv([1, -1], conv([1, -p2], [1, -p3]))
%% Matched Z transform
% directly map poles and zeros using z<= e^(-sT)

T = pi/10;

% Continuous-time poles
p1 = 1;
p2 = exp((-1-1j)*T);
p3 = exp((-1+1j)*T);

disp(['p1 = ', num2str(p1), ', p2 = ', num2str(p2), ', p3 = ', num2str(p3)]);

% Create the denominator polynomial
denom = conv([1, -1], conv([1, -p2], [1, -p3]));

disp(['denom = ', num2str(denom)]);

% Create the discrete-time transfer function
gol_d2 = tf([1, 0, 0, 0], real(denom), T);

figure;
pzmap(gol_d2);

% Calculate the gain factor
impulse_response_gol = impulse(gol, [0:100]*T);
impulse_response_gol_d2 = impulse(gol_d2, [0:100]*T);
gain = impulse_response_gol(end) / impulse_response_gol_d2(end);
disp(['Gain = ', num2str(gain)]);

%% Compare Impulse Responses
time_axis = (0:29)*T;
figure()
[yout_c1, time_out1] = impulse(gol, time_axis);
hold on
[yout_d2, time_out2] = impulse(gol_d, time_axis);
hold on
[yout_d3, time_out3] = impulse(gol_d2, time_axis);

figure;
title("Comparative Impulse Responses");
plot(time_out1, yout_c1, 'LineWidth', 1.5, 'DisplayName', 'Continuous-Time');
hold on;
plot(time_out2, yout_d2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Impulse Invariance');
plot(time_out3, gain*yout_d3, 'LineWidth', 1.5, 'DisplayName', 'Matched-Z');
hold off;
legend;
grid;

figure;
bode(gol);
hold on
bode(gol_d);
hold on
bode(gol_d2);
legend('Continuous-Time', 'Impulse Invariance', 'Matched-Z');
title('Bode Plot');



%% Power control Loop
K = 1;
golz = tf(0.83155 * K, [1, -1, 0], 1);
disp(golz);

figure;
nyquist(golz);
axis([-1.4, 0.5, -10, 10]);

figure;
nyquist(-1/golz);
axis([-2.5, 1.5, -2.5, 2.5]);
axis equal;

[rlist, klist] = rlocus(golz);
figure;
plot(real(rlist), imag(rlist));
title('Root Locus');
angle = linspace(0, 2*pi, 512);
hold on;
plot(real(exp(1j*angle)), imag(exp(1j*angle)));
hold off;
axis([-2, 2, -1.5, 1.5]);


sys1 = tf(0.01464, [1, -1], 1);
sys2 = tf(56.8, [1, 0], 1);
n = 0:24;

% Step is done on the closed-loop system
K = 0.5;
figure;
[yout, n] = step(1500*feedback(K*sys1, sys2), n);
stairs(n, squeeze(yout), 'DisplayName', sprintf('K = %.2f', K));

K = 0.25;
[yout, n] = step(1500*feedback(K*sys1, sys2), n);
hold on;
stairs(n, squeeze(yout), 'r', 'DisplayName', sprintf('K = %.2f', K));
hold off;

grid on;
legend;
title('Discrete Power Control Loop Step Response');
ylabel('Power Gain [dB]');
xlabel('Sample Number');