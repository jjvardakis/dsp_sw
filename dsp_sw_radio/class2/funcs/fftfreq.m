function f = fftfreq(N, d)
    % FFTFREQ Return the FFT sample frequencies.
    %
    %   f = FFTFREQ(N) returns a vector of length N containing the frequency
    %   values of the FFT output. The returned frequency values are in cycles per
    %   unit of the sampling interval.
    %
    %   f = FFTFREQ(N, d) specifies the sampling interval (inverse of the
    %   sampling frequency). If not specified, d defaults to 1.

    if nargin < 2
        d = 1;
    end

    if mod(N, 2) == 0
        % For even N
        f = (-N/2:N/2-1) / (N * d);
    else
        % For odd N
        f = (-(N-1)/2:(N-1)/2) / (N * d);
    end
end