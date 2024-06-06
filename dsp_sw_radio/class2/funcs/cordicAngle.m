function [I1, Q1] = cordicAngle(I, Q, angle, N)
    % Rotates vector I, Q by angle.
    % I, Q: Initial vector real and imaginary values,
    % angle: Angle in radians within +/- pi/2,
    % N: Number of iterations.
    % Returns rotated I, Q.
    binaryWord = binAngle(angle, N);
    [I1, Q1] = cordicSequence(I, Q, binaryWord);
end