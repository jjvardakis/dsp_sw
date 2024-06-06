function result = pulseShape(symbols, oversamp, alpha, mode, sym_length)
    % returns a pulse shaped waveform defined by list of complex symbols with output
    % time aligned with input (zero-phase pulse shape filtering), rather than predict the output
    % of an FIR pulse shaping filter implementation which would be delayed due to the filter.
    %
    % oversamp is the oversampling factor (integer >=2)
    % alpha is the filter roll-off factor (between 0 < alpha <= 1)
    % type: 'rrc' is root-raised cosine, 'rc' is raised cosine
    % sym_length: length of the filter impulse response in number of symbols (default = 46) 
    
    % (Filter is linear phase with delay = half the length of the filter)
    % prepend symbols with first symbol by half length of filter to pre-fill
    % and append with last symbol by half length to process all symbols
    sym_prepend = [symbols(1) * ones(1, sym_length/2), symbols];
    syms = [sym_prepend, symbols(end) * ones(1, sym_length/2)];
    % zero-insert upsampling
    tx = zeros(1, length(syms) * oversamp);
    tx(1:oversamp:end) = syms*oversamp;
    % apply pulse shape filter
    result = filter(psCoeff(oversamp, alpha, mode, sym_length), 1, tx);
    
    % truncate first samples due to prepend and append to align output with input
    result = result(sym_length * oversamp + 1:end);
end