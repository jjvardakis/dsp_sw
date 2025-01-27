clear all
close all
clc

fund = sin(2*pi*.01*[1:100]);
second = .5*sin(4*pi*.01*[1:100]);
third = .3*sin(6*pi*.01*[1:100]);

figure()
subplot(2,1,1)
plot(fund)
hold on
plot(second)
hold on
plot(fund+second)
legend("first", "second", "Sum")
xlabel("Time (samples)")
title("First and Second Harmonic - Assymetric")

subplot(2,1,2)
plot(fund)
hold on
plot(third)
hold on
plot(fund+third)
legend("first", "third", "Sum")
xlabel("Time (samples)")
title("First and Third Harmonic - Symmetric")
% tight_layout()


% Parameters
M = 16; % Number of symbols in the constellation
BW = 1e3; % rate at which each symbol is sent


osf1 = 8; % samples per symbol
fs = BW * osf1;

alpha = 0.1;   %  Roll-off factor (set between 0 and 1)
rc_span = 36; % Duration of impulse response
h_rrc = rcosdesign(alpha, rc_span, osf1,"sqrt")/sqrt(osf1);
t = linspace(-rc_span/2,rc_span/2, rc_span * osf1+1);

h_rc =conv(h_rrc,h_rrc)/osf1;
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
bits = 12;
agc  = -10;   % agc level in dBFS DC
osf1 = 8; % oversampling factor of 8

% create dictionary for mapping data to symbols
const_map = qammod(data,M);
tx = upsample(const_map,osf1);
h_rc = conv(h_rrc,h_rrc)/osf1;
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

% consider an SDR signal sampled with sampling rate >> bw

% we can emulate this by resampling the tx_shaped waveform above to a higher rate and then selecting 
% a higher image of the interpolated images to have an upsampled and digital IF waveform:


I = 10;   % interpolation factor
fs_x10 = fs*I;

% zero insert
txhr = zeros(length(tx_shapedquant1) * I,1);
txhr(1:I:end) = tx_shapedquant1;   %also increased signal level relative to reference used
[txhr_fd, faxis] = win_fft(txhr.',fs*I);
faxis_centered = faxis-fs_x10/2;
%txhr = sig.resample(tx_shaped, len(tx_shaped) * 10)

figure()
N = length(txhr);
plot(faxis_centered,20*log10(abs(fftshift(txhr_fd))))

% create multiband filter targets:
bandcenters = [1:4]*2*fs/(fs*I);
width = .9/(fs/2*I)*10^3;
bandedges = [0, width];
for center=bandcenters
    bandedges=[bandedges center-width center+width];
end
% Append two more values to represent the last band
bandedges = [bandedges, 5*2*fs /(fs*I) - width, 5*2*fs /(fs*I)];

% create complex filter to select single image by heterodying low pass coeff of multiband filter:
ntaps =91;

% Generate the filter using firls
band_select = firls(ntaps-1, bandedges, [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

coeff = band_select .* exp(2j*pi/I * [1:ntaps]);
h = freqz(coeff,1,2^16, 'whole');
digif = real(filter(coeff,1,txhr));

sigma = 1e-5;
noise = sigma * (randn(length(digif),1)+1i*randn(length(digif),1));

digifsig = digif + noise;

figure()
[fout, faxis] = win_fft(digifsig.',fs_x10);


figure
plot(linspace(-fs/2*I+1,fs/2*I,2^16), db(abs(fftshift(h))));   % superimposing filter response
hold on
plot(faxis_centered, db(abs(fftshift(txhr_fd))))
hold on
plot(faxis-fs*I/2,20*log10(abs(fftshift(fout))))
legend("Baseband signal","Multi-band filter","IF+noise signal")
title("Selection of 16-QAM Digital IF Spectrum");

%% Example 4: Non-linear transfer function of the magnitude %%

% demonstrate
vin = linspace(-5,5,1000);
vout = compress(vin, 2);
vout2 = compress(vin, 3, false);
figure()
plot(vin,vout)
hold on
plot(vin,vout2)
title("Vin vs Vout Showing Non-linear Compression")
xlabel("Vin")
ylabel("Vout")
grid on

third_sat = 0.25;
second_sat = 0.1;

% Adds 3rd order harmonics
sigout = compress(digifsig, third_sat);

% Adds 2nd order harmonics
sigout = compress(sigout, second_sat, false, true);

figure;
[Pxx, f] = win_fft(sigout.', fs_x10);
plot(f-fs_x10/2, 20*log10(abs(fftshift(Pxx))), 'LineWidth', 2);
title('16-QAM Digital IF Spectrum w/ Nonlinear Distortion');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

figure()
vin = linspace(-.2,.2,100);
vout = distort(vin,third_sat,second_sat);
figure()
plot(vin,vout)
title("Actual Non-Linear Distortion Used")
xlabel("Vin")
ylabel("Vout")
grid on

t = (0:length(sigout)-1) * 1/(fs_x10);

tone = std(digifsig)* sqrt(2) * cos(2*pi*4e3*t) + noise.';

% % Adds 3rd order harmonics
% tonedistort = compress(tone.', third_sat);
% 
% % Adds 2nd order harmonics
% tonedistort = compress(tonedistort.', second_sat, false, true);
tonedistort = distort(tone.',third_sat,second_sat);

figure;
[Pxx, f ] = win_fft(tonedistort.', fs_x10);
plot(f-fs_x10/2, 20*log10(abs(fftshift(Pxx))), 'LineWidth', 2);
title('Single Tone at Same Power with the Same Nonlinear Distortion');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

%%%%%%%%%%%%%
t = (0:length(sigout)-1) * 1/(fs_x10);

tone1 = std(digifsig) * cos(2*pi*3.8e3*t) + noise.';
tone2 = std(digifsig) * cos(2*pi*4.2e3*t) + noise.';
tonedistort2 = distort((tone1 + tone2).',third_sat,second_sat);

figure;
[Pxx, f ] = win_fft(tonedistort2.', fs_x10);
plot(f-fs_x10/2, 20*log10(abs(fftshift(Pxx))), 'LineWidth', 2);
title('Two Tones at Same Power with the Same Nonlinear Distortion');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

% Define the distortion function
function distorted = distort(waveform,third_sat,second_sat)
    % Same operation as above adding even and odd order distortion
    distorted = compress(compress(waveform, third_sat), second_sat, false, true);
end
function db_values = db(x)
    % Replace zeros with the smallest positive normal value
    x_safe = x;
    x_safe(x == 0) = realmin('double');
    
    % Calculate dB values
    db_values = 20 * log10(abs(x_safe));
end

function compressed = compress(vin, limit, pos, neg)
% will just add a hyperbolic tangent compression while keeping same gain for small
% signals and hardlimit at limit, compression to both positive and negative
% values when pos=True and neg=True (default)
% Default values
if nargin < 4
    neg = true;
end
if nargin < 3
    pos = true;
end
if nargin < 2
    limit = 1;
end
if ~(pos || neg)
    % doesn't do anything
    compressed = vin;
elseif ~pos
    compressed = vin .* (vin < 0) + limit * tanh(vin/limit) .* (vin >= 0);
elseif ~neg
    compressed = vin .* (vin >= 0) + limit * tanh(vin/limit) .* (vin < 0);
else
    compressed = limit * tanh(vin/limit);
end
end