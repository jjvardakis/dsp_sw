function result = db(x)
    % returns dB of number and avoids divide by 0 warnings
    x_safe = max(x, eps(class(x))); % avoids divide by 0
    result = 20 * log10(abs(x_safe));
end