function poly = spPoly(num, den, x, prec)
    % return a symbolic polynomial or polynomial fraction
    % x is string of dependent variable (defaults to 's')
    % num is array-like coefficients of numerator polynomial
    % den (optional) is array-like coefficients of denominator polynomial
    % prec is accuracy of the decimal digits used (default 6)
    
    if nargin < 4
        prec = 6;
    end

    syms x;

    numRound = round(num, prec);
    denRound = round(den, prec);
    
    poly = poly2sym(numRound, x) / poly2sym(denRound, x);
end
