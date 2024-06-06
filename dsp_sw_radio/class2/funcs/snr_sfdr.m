function [snr, sfdr] = snr_sfdr(waveform)
    N = length(waveform);
    
    win = kaiser(N, 16);
    
    % Remove DC offset
    alpha = 0.99;
    waveform_nodc = filter([(1+alpha)/2, -(1+alpha)/2], [1, -alpha], waveform);
    
    % Main lobe with beta = 16 ~ 10.4: use +/-6 bins from peak as "signal"
    fout = fft(waveform_nodc .* win);
    
    % Signal selection
    [~, ksig] = max(abs(fout));
    sig_indices = mod(ksig + (-6:6), N) + 1;

    % Remove beyond signal for the purpose of measuring noise
    excise_indices = mod(ksig + (-20:20), N) + 1;
    
    % Sum square signal components
    signal = sum(abs(fout(sig_indices)).^2);
    
    % Remove signal
    no_sig = fout;
    no_sig(excise_indices) = 0;
  
    % Sum noise components
    noise = sum(abs(no_sig).^2);
    
    % Find the strongest spur
    [~, kspur] = max(abs(no_sig));
    spur_indices = mod(kspur + (-6:6), N) + 1;
    
    % Sum square spur components
    spur = sum(abs(no_sig(spur_indices)).^2);
    
    snr = 10*log10(abs(signal/noise));
    sfdr = 10*log10(abs(signal/spur));
end