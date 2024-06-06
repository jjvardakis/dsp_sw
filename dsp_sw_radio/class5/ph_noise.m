function [F, P, pnoise] = ph_noise(signal, Fs, fpoints, npowers)
% MATLAB equivalent of ph_noise function

    % Compute power spectral density of the signal
    [Pxx, F] = pwelch(signal, [], [], [], Fs);

    % Initialize the phase noise
    pnoise = zeros(size(signal));

    % Add phase noise at specified frequencies
    for i = 1:length(fpoints)
        % Find the index corresponding to the frequency point
        [~, idx] = min(abs(F - fpoints(i)));

        % Add phase noise
        phi = sqrt(10^(npowers(i) / 10) * (F(idx + 1) - F(idx)));
        pnoise = pnoise + phi * randn(size(signal));
    end

    % Apply the phase noise to the signal
    signal_with_noise = signal .* exp(1j * pnoise);

    % Compute the power spectral density of the signal with noise
    [P, F] = pwelch(signal_with_noise, [], [], [], Fs);
end