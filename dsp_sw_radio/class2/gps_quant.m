% Load captured data as done in Notebook 1, although
% here we assign each data sample to be +1 or -1 (hard limit or 1 bit A/D converter)

% constants

f_s = 20e6;        % sampling rate
f_chip = 1.023e6;  % chip rate
freq_if = 2.5e6;   % digital IF frequency
doppler = 1.088e3; % net frequency offset due to Doppler and clock offsets in Hz for SV24
phase = 2.985;    % initial phase offset for SV24 in radians
spc = f_s/f_chip;  % number of samples per chip

% load captured signal
% GPS Signal at 1575.42 MHz downconverted to the digital IF frequency
% and sampled (using Lecroy scope) at 20 MSPS 8 bits
% on April 7, 2007 at 1:35pm in Woburn MA
load("E:\DSP_course_for_Software_Radios\scripts\python scripts\class1\data\gps\gps_data.txt");
data = size(gps_data);
for i = 1:length(gps_data)
    if(gps_data(i)>0)
        data(i) = 1;
    else
        data(i) = -1;
    end
end

n = length(data);

% create time vector
t = [1:length(data)]/f_s;


figure()
plot(t,data)
xlabel('Time [s]')
ylabel('Magnitude [v]')
title(sprintf("Captured Data, %d samples",n))


% downconvert signal from digital IF + Doppler to baseband:
f_shift = 2*pi*(-freq_if - doppler);
baseband= data .* exp(1i * (f_shift*t - phase));

N = 1023;
lag = [-N+1:N];

% create C/A code for SV24
prn24 = cacode(24,1);    % will generate 1023 chips (1 code sequence)

% confirm 1023 chips: typical way to count % of items in generator:
n_chips = length(prn24);
fprintf("Number of chips is %d\n",n_chips)
for i = 1:length(prn24)
    prn24(i) = prn24(i)*2-1;    % will generate 1023 chips (1 code sequence) mapped to +1 / -1
end

% resample to fs
prn24_fs = resample(prn24, round(f_s/n_chips), round(f_chip/n_chips));

b = fft(baseband);
a = fft(prn24_fs,length(b));
c = ifft(a .* conj(b));

figure()
plot(real(c))
title("Correlated Symbols for SV24")
xlabel('Time [s]')
ylabel('Magnitude')