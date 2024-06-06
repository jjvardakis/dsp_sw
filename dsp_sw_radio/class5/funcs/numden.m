function [num, den] = numden(tf)
    [num, den] = tfdata(tf);
    num = cell2mat(num);
    den = cell2mat(den);
end
