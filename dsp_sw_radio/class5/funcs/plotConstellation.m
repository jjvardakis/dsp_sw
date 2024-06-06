function samps = plotConstellation(tx, offset, title, oversamp, interpol)
    % Resample to higher rate for purpose of plotting
    tx_upsamp = interp(tx, interpol);

    % Omit first 30 symbols and last 100 samples
    tx_upsamp = tx_upsamp(oversamp * interpol * 30:end - 100);

    figure;
    plot(real(tx_upsamp), imag(tx_upsamp));

    % Offset to select the correct timing sample and skip the first 30 symbols
    % (waveform is now sampled 40 samples/symbol for the purpose of plotting)
    samps = real(tx_upsamp(offset:oversamp * interpol:end)) + 1j * imag(tx_upsamp(offset:oversamp * interpol:end));

    hold on;
    plot(real(samps), imag(samps), 'r.');
    hold off;
    axis equal;
end