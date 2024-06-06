function [I1, Q1] = cordicSequence(I, Q, sequence)
    % Rotates vector I, Q using +/- directions given in sequence.
    % I, Q: Initial vector real and imaginary values,
    % sequence: Iteration direction as a binary string (e.g., '10110' for a 5 iteration CORDIC).
    % Returns rotated I, Q.
    I1 = I;
    Q1 = Q;
    for n = 1:length(sequence)
        d = 2*str2double(sequence(n))-1;
        I2 = I1 - d*Q1/(2^n);
        Q2 = d*I1/(2^n) + Q1;
        I1 = I2;
        Q1 = Q2;
    end
end