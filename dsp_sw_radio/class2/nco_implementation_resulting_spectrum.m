
%% Example 1 
fout = 1.132e6;
fclk = 10e6;
acc_s = 24;

fres = fclk / 2^acc_s;

fcw = round(fout / fres);

fprintf('FCW: %d\n', fcw);

% Initialize the NCO generator
nco_gen = nco(fcw, 0, acc_s, 14, 12, 2^15, 0);

% NCO generator is a generator, so use a loop to get the results
result = nco_gen;

% frequency spectrum of NCO

figure
plot(result(1:50))

figure()
% using cusomized fft module imported earlier
[Pxx,f] = win_fft(result/max(result), fclk)
plot(f-fclk/2, 20*log10(abs(fftshift(Pxx))), 'LineWidth', 2);
title('NCO Output Spectrum')


%% Example 2: NCO Dithering example
nco_gen_dither = nco(fcw, 0, acc_s, 14, 12, 2^15, 6);

result_dither = nco_gen_dither;

figure()
% using cusomized fft module imported earlier
[Pxx,f] = win_fft(result_dither/max(result_dither), fclk)
plot(f-fclk/2, 20*log10(abs(fftshift(Pxx))), 'LineWidth', 2);
title('NCO Output Spectrum - Dither Enabled')

% Frequency sweep with NCO

fout = linspace(1e3, 15e3, 2^14);
fclk = 10e6;
acc_s = 32;


fres = fclk / 2^acc_s;

fcw = round(fout / fres)

nco_gen = nco(fcw, 0, acc_s, 14, 12, 2^15, 0);

figure()
plot(nco_gen)

 
%% Example 3: NCO Spurs and SNR
N = 2^12;
wink = kaiser(N, 16);
winr = ones(1, N);

figure;

pad = 2^18;
faxis = fftfreq(pad);

plot(faxis, db(fftshift(fft(winr, pad))/N), 'LineWidth', 1.5, 'DisplayName', 'Rectangular');
hold on
plot(faxis, db(fftshift(fft(wink, pad))/N), 'LineWidth', 1.5, 'DisplayName', 'Kaiser');
axis([-0.005, 0.005, -200, 0]);
grid on;
xlabel('Normalized Frequency (rad/sample)');
ylabel('dB');
title('Comparing Kernel of Rectangular Window to Kaiser Window');
legend;


%% Example 4: Test case for NCO
acc_s = 16;
lut_in = 16;
lut_out = 20;
fclk = 10e6;
fstep = fclk / 2^acc_s;
fprintf('Frequency step size is %.2e Hz\n', fstep);

nsamps = 2^15;

% Create NCO generator
nco_gen = Cnco(acc_s, lut_in, lut_out);

% Case 1: Integer Submultiple of the sampling rate
fcw1 = 2^12;
fcw = ones(1, nsamps) * fcw1;
pcw = zeros(1, nsamps);

% Prime NCO generator
result1 = zeros(1, nsamps);
for i = 1:nsamps
    result1(i) = nco_gen(fcw(i), pcw(i),0);
end

% Filter with DC blocking filter
alpha = 0.99;
result1 = filter([(1+alpha)/2, -(1+alpha)/2], [1,-alpha], result1);

% Case 2: Non-Integer Submultiple of the sampling rate
fcw2 = 2^12 + 37;
fcw = ones(1, nsamps) * fcw2;
pcw = zeros(1, nsamps);

result2 = zeros(1, nsamps);
for i = 1:nsamps
    result2(i) = nco_gen(fcw(i), pcw(i), 0);
end

% Filter with DC blocking filter
result2 = filter([(1+alpha)/2, -(1+alpha)/2], [1,-alpha], result2);

% Plot the results
figure('Position', [100, 100, 800, 300]);

subplot(1, 2, 1);
[Pxx,f] = win_fft(result1 / max(result1), fclk);
plot(f-fclk/2, 20*log10(abs(fftshift(Pxx))), 'LineWidth', 2);
title(['FCW = $2^9$']);

subplot(1, 2, 2);
[Pxx,f] = win_fft(result2 / max(result2), fclk);
plot(f-fclk/2, 20*log10(abs(fftshift(Pxx))), 'LineWidth', 2);
title(['FCW = ', num2str(fcw2)]);

sgtitle('NCO Test Case');

% Calculate SNR and SFDR for the first test case
[y, sfdr] = snr_sfdr(result1.');
fprintf('SNR for FCW = %d is %.2f dB and SFDR is %.2f dB\n', fcw1, y, sfdr);

% Calculate SNR and SFDR for the second test case
[y, sfdr] = snr_sfdr(result2.');
fprintf('SNR for FCW = %d is %.2f dB and SFDR is %.2f dB\n', fcw2, y, sfdr);


%% DC nulling filter
alpha = 0.99;
% [b, a] = fir1(64, 1, 'DC-0');
a = [(1+alpha)/2, -(1+alpha)/2];
b = [1, -alpha];
[h, w] = freqz( a, b, 2^16, "whole");

figure;
plot(w, db(h));
% axis([0, pi, -100, 10]);
grid on;
xlabel('Normalized Frequency');
ylabel('dB');
title('DC Nulling Filter');

N = length(result2);
win = kaiser(N, 16);

fout = fft(result2.' .* win);

ksig = find(fout == max(fout), 1);
% main lobe with beta = 16 ~ 10.4: use +/-6 bins from peak as "signal"
sig_indices = mod(ksig + (-6:6), N) + 1;

% remove beyond signal for the purpose of measuring noise
excise_indices = mod(ksig + (-20:20), N) + 1;

% sum square signal components
signal = sum(abs(fout(sig_indices)).^2);

% remove signal    
no_sig = fout;
no_sig(excise_indices) = 0;
kspur = find(no_sig == max(no_sig), 1);
spur_indices = mod(kspur + (-6:6), N) + 1;

figure;
plot(db(fout));
hold on;
plot(db(no_sig));
axis([-500, N+500, -20, 200]);
title("Demo of Signal Extraction for SNR measurement");
grid on;
ylabel('dB');
xlabel("Frequency bin");
legend('Original Spectrum', 'Spectrum with Signal Removed');


%% sweep through all FCW's and report SNR and SFDR: (this takes a while to run for large FCW ranges!)

acc_s = 16;
nco_gen = Cnco(acc_s, 14, 12);

snrs_sfdrs = [];
fcw_start = 200;
fcw_stop = 500;
heart_beat = 50;
fprintf('Sweeping from FCW = %d to FCW = %d\n', fcw_start, fcw_stop);
fprintf('Iteration number (prints every %d iters):\n', heart_beat);

fcw_values = 0:(2^acc_s - 1);

for fcw_index = fcw_start:fcw_stop
    fcw_value = fcw_values(fcw_index + 1); % MATLAB indexing starts from 1
    fcw = ones(1, 2^15) * fcw_value;
    pcw = zeros(1, 2^15);
    
    if mod(fcw_value, heart_beat) == 0
        fprintf('%d\n', fcw_value);
    end
    
    result = zeros(1, 2^15);
    for i = 1:2^15
        [freq, phase] = deal(fcw(i), pcw(i));
        result(i) = nco_gen(freq, phase,0);
    end
    snr_sfdr(result.');
    [y, sfdr] = snr_sfdr(result.');
    snrs_sfdrs= [snrs_sfdrs; [y, sfdr]];
end

snrs = snrs_sfdrs(:, 1);
sfdrs = snrs_sfdrs(:, 2);

figure;
plot(fcw_values(fcw_start:fcw_stop), sfdrs, 'LineWidth', 2, 'DisplayName', 'SFDR');
hold on;
plot(fcw_values(fcw_start:fcw_stop), snrs, 'LineWidth', 2, 'DisplayName', 'SNR');
grid on;
xlabel('FCW');
ylabel('dB');
title(sprintf('SNR and SFDR vs FCW for %d-bit Accum, %d-bit LUT in, %d-bit LUT out', acc_s, lut_in, lut_out));
legend('Location', 'Best');