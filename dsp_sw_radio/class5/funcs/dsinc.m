function result = dsinc(x)
    % derivative sinc function (careful to scale properly as dx is not included other than pi if it also has a derivative)
    num = pi * x .* cos(pi * x) - sin(pi * x);
    den = pi * x.^2;
    result = num ./ den;
    
    % Handle division by zero
    result(isnan(result) | isinf(result)) = 0;
end