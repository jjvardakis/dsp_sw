% 16QAM constellation already created above
num_symb = 2^12;

data = floor(16*rand(num_symb,1));

% create dictionary for mapping data to symbols
const_map = qammod(data,16);
% Normalize power of constellation
qam16_const = const_map / std(const_map);

qam16_peak = max(abs(qam16_const));

% 16APSK constellation: 
% (example parameters from DBS-S2: 
% https://www.etsi.org/deliver/etsi_en/302300_302399/30230701/01.04.01_20/en_30230701v010401a.pdf
gamma = 2.75;
apsk16 = apskmod(data,[4 12],[1,gamma], [pi/4 0]); 
% Normalize power of constellation
apsk16_const = apsk16 / std(apsk16);

apsk16_peak = max(abs(apsk16_const));

figure;
subplot(1,2,1);
plot(real(qam16_const), imag(qam16_const), 'r.');

hold on;
plot(qam16_peak * cos(linspace(0,2*pi,100)), qam16_peak * sin(linspace(0,2*pi,100)), 'b--', 'LineWidth', 0.3);
hold off;

title('16QAM');
axis([-1.5, 1.5, -1.5, 1.5]);
subplot(1,2,2);
plot(real(apsk16_const), imag(apsk16_const), 'r.');

hold on;
plot(apsk16_peak * cos(linspace(0,2*pi,100)), apsk16_peak * sin(linspace(0,2*pi,100)), 'b--', 'LineWidth', 0.3);
hold off;

axis([-1.5, 1.5, -1.5, 1.5]);
title('16APSK');

fprintf('Peak vector for 16QAM = %.3f\n', qam16_peak);
fprintf('Peak vector for 16APSK = %.3f\n', apsk16_peak);
fprintf('Min distance for 16QAM = %.3f\n', abs(qam16_const(1) - qam16_const(2)));
fprintf('Min distance for 16APSK (inner ring) = %.3f\n', abs(apsk16_const(1) - apsk16_const(2)));
fprintf('Min distance for 16APSK (outer ring) = %.3f\n', abs(apsk16_const(5) - apsk16_const(6)));
