function result = psCoeff(oversamp, alpha, mode, sym_length)
    % returns coefficients for FIR root-raised-cosine or raised-cosine pulse-shaping filter
    if strcmp(mode, 'rrc')
        pulse = rcosdesign(alpha, sym_length, oversamp, 'sqrt');
    else
        pulse = rcosdesign(alpha, sym_length, oversamp, 'normal');
    end
    
    index = linspace(-sym_length/2, sym_length/2, sym_length * oversamp);
    
    % scaling to normalize impulse response for use as filter with 0 dB passband
    result = pulse / oversamp;
end