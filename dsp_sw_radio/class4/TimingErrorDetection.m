% Eye Diagram

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
