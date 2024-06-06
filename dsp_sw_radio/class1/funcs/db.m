function db_values = db(x)
    % Replace zeros with the smallest positive normal value
    x_safe = x;
    x_safe(x == 0) = realmin('double');
    
    % Calculate dB values
    db_values = 20 * log10(abs(x_safe));
end