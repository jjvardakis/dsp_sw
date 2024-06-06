function output = nco(fcw, pcw, acc_size, lut_in_size, lut_out_size, nsamp, dither)
    acc = pcw;

    % Convert input_values to an iterator if a scalar or list-like input
    if ~isa(fcw, 'cell') && numel(fcw) == 1
        % Input is a scalar
        fcw = repmat(fcw, 1, nsamp);
    end

    output = zeros(size(fcw));
    
    for count = 1:numel(fcw)
        if nsamp && count > nsamp
            break;
        end
        samp = fcw(count);
        acc = mod(samp + acc + round(2^(dither-1) * randn(1)), 2^acc_size);
        lut_in = floor(acc / 2^(acc_size-lut_in_size));
        angle = lut_in / 2^lut_in_size * 2 * pi;
        output(count) = round(2^(lut_out_size-1) * cos(angle));
    end
end