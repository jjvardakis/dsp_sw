%% Quantization effects
% Parameters
M = 16; % Number of symbols in the constellation
BW = 1e3; % rate at which each symbol is sent


osf1 = 8; % samples per symbol
fs = BW * osf1;

alpha = 0.1;   %  Roll-off factor (set between 0 and 1)
rc_span = 36; % Duration of impulse response
h_rrc = rcosdesign(alpha, rc_span, osf1,"sqrt")/sqrt(osf1);
t = linspace(-rc_span/2,rc_span/2, rc_span * osf1+1);
disp(length(h_rrc))

h_rc =conv(h_rrc,h_rrc);
[h, w] =freqz(h_rc, 1, "half",fs)

% Plot frequency response
figure;
plot(w, 20*log10(abs(h)));
title('Frequency Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;


% create random data source
num_symb = 2^12;

data = floor(16*rand(num_symb,1));
bits = 12
agc  = -10   % agc level in dBFS DC

% create dictionary for mapping data to symbols
const_map = qammod(data,M);
tx = upsample(const_map,osf1);
h_rc = conv(h_rrc,h_rrc);
tx_shaped = filter(h_rc,1,tx);

tx_shaped = tx_shaped/std(tx_shaped) *10^(agc/20);
tx_shapedi = round(2^(bits-1)*real(tx_shaped))/2^(bits-1);
tx_shapedq = round(2^(bits-1)*imag(tx_shaped))/2^(bits-1);

tx_shapedquant1 = complex(tx_shapedi, tx_shapedq);

% compute dBFS for rms level:
dBFStx = db(std(tx_shapedquant1));
sprintf("RMS signal is %.2f dBFS DC", dBFStx)

% compute peak signal for this dataset
sprintf("Peak signal is %.1fdb dBFS DC",db(max(abs(tx_shapedquant1))))

% compute quantization noise
error = tx_shapedquant1-tx_shaped;
sprintf("Total Quantization Noise is %.1f dBfS DC", 20*log10(std(error)))

%% plot the frequency spectrum of the signal along with its papr
figure()
N = length(tx_shapedquant1);
plot(fs/N*(-N/2:N/2-1),20*log10(abs(fftshift(fft(tx_shapedquant1)))))
title(sprintf("%.1f KHz 16-QAM Spectrum, fs=%.1f KHz,  %d bits", BW/1000, fs/1000, bits))
% xlim([-fs*oversamp1/2000 fs*oversamp1/2000])

[papr, results] = ccdf(tx_shaped);

figure('Position', [100, 100, 600, 400]);
semilogy(papr, results);
xlabel('PAPR (dB)');
ylabel('% exceeding PAPR');

% Use sprintf to create a formatted string for the title
title_str = sprintf('CCDF for Single Carrier 16-QAM, RRC \\alpha = %.2f', alpha);
title(title_str);
grid on;

%% additional upsampling 
oversamp2 = 5;
oversamp = osf1*oversamp2;
fs = fs * oversamp2;
tx = upsample(tx_shapedquant1,oversamp2);
h_rc = conv(h_rrc,h_rrc);
tx_shaped = filter(h_rc,1,tx);

tx_shaped = tx_shaped/std(tx_shaped) *10^(agc/20);
tx_shapedi = round(2^(bits-1)*real(tx_shaped))/2^(bits-1);
tx_shapedq = round(2^(bits-1)*imag(tx_shaped))/2^(bits-1);

tx_shapedquant2 = complex(tx_shapedi, tx_shapedq);

% compute dBFS for rms level:
dBFStx = db(std(tx_shapedquant2));
sprintf("RMS signal is %.2f dBFS DC", dBFStx)

% compute peak signal for this dataset
sprintf("Peak signal is %.1fdb dBFS DC",db(max(abs(tx_shapedquant2))))

% compute quantization noise
error = tx_shapedquant2-tx_shaped;
sprintf("Total Quantization Noise is %.1f dBfS DC", 20*log10(std(error)))

%% plot the frequency spectrum of the signal along with its papr
figure()
N = length(tx_shapedquant2);
plot(fs/N*(-N/2:N/2-1),20*log10(abs(fftshift(fft(tx_shapedquant2)))))
title(sprintf("%.1f KHz 16-QAM Spectrum, fs=%.1f KHz,  %d bits", BW/1000, fs/1000, bits))
% xlim([-fs*oversamp/2000 fs*oversamp/2000])
[papr, results] = ccdf(tx_shapedquant2);

figure('Position', [100, 100, 600, 400]);
semilogy(papr, results);
xlabel('PAPR (dB)');
ylabel('% exceeding PAPR');

% Use sprintf to create a formatted string for the title
title_str = sprintf('CCDF for Single Carrier 16-QAM, RRC \\alpha = %.2f', alpha);
title(title_str);
grid on;


figure 
histogram(real(tx_shapedquant2),500)
title("Histogram for real axis of tx waveform")

eyediagram(tx_shapedquant1(1:2048),2*osf1)

%% 
tx_resamp = resample(tx_shapedquant1,8,1);
figure()
plot(real(tx_resamp(581:1800)))

interp =16;
upsampled = resample(tx_shapedquant1,16,1);
fs = fs*interp;
offset = 8;
tx_offset = upsampled(offset:interp:end).*exp(1i*.3);
eyediagram(tx_offset(1:2048),2*osf1)

figure()
scatter(real(tx_offset),imag(tx_offset))
hold on
start = 2*osf1+1;
scatter(real(tx_offset(start:8:end-start)),imag(tx_offset(start:8:end-start)))


back_off = -12              % backoff from full scale in dB
n_bits = 6                  % number of bits for I and for Q

% normalized to 1 = dBFS DC
tx_offset = tx_offset/std(tx_offset) * 10^(back_off/20);

tx_shapedi = round(2^(n_bits-1)*real(tx_offset))/2^(n_bits-1);
tx_shapedq = round(2^(n_bits-1)*imag(tx_offset))/2^(n_bits-1);

tx_offset_quant = tx_shapedi + 1j*tx_shapedq;

figure()
scatter(real(tx_offset_quant), imag(tx_offset_quant))
hold on
scatter(real(tx_offset_quant(121:osf1:end)),imag(tx_offset_quant(121:osf1:end)))

figure()
N = length(tx_offset_quant);
plot(fs/N*(-N/2:N/2-1),20*log10(abs(fftshift(fft(tx_offset_quant)))))
title(sprintf("%.1f KHz 16-QAM Spectrum, fs=%.1f KHz,  %d bits", BW/1000, fs/1000, bits))

% compute dBFS for rms level:
dBFStx = db(std(tx_offset_quant));
sprintf("RMS signal is %.2f dBFS DC", dBFStx)

% compute peak signal for this dataset
sprintf("Peak signal is %.1fdb dBFS DC",db(max(abs(tx_offset_quant))))


function db_values = db(x)
    % Replace zeros with the smallest positive normal value
    x_safe = x;
    x_safe(x == 0) = realmin('double');
    
    % Calculate dB values
    db_values = 20 * log10(abs(x_safe));
end

function [papr, ccdf_values] = ccdf(x)
    papr = 0:0.1:12;  % dB vector of PAPR values
    power = abs(x).^2;
    power_ratio = power / mean(power);
    pdb = 10 * log10(power_ratio);
    
    ccdf_values = zeros(size(papr));
    for i = 1:length(papr)
        ccdf_values(i) = sum(pdb > papr(i)) / length(x);
    end
end