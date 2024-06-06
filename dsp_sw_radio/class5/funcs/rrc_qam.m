function [result,symbols] = rrc_qam(num_symb, order, oversamp, alpha, sym_length)
    % sym_length: duration of impulse response
    % oversamp: samples per symbol
    % alpha: Roll-off factor (set between 0 and 1)

    % Set default values if not provided
    if nargin < 5
        sym_length = 46;
    end
    
    if nargin < 4
        alpha = 0.25;
    end
    
    if nargin < 3
        oversamp = 4;
    end
    
    symbols = qam(num_symb, order);
    result = pulseShape(symbols, oversamp, alpha, 'rrc', sym_length);
end