function result = rc_qam(num_symb, order, oversamp, alpha, sym_length)
    % sym_length: duration of impulse response
    % oversamp: samples per symbol
    % alpha: Roll-off factor (set between 0 and 1)
    symbols = qam(num_symb, order);
    result = pulseShape(symbols, oversamp, alpha, 'rc', sym_length);
end