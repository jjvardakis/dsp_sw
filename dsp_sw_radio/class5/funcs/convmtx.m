function result = convmtx(h, n)
    % creates the convolution (Toeplitz) matrix, which transforms convolution to matrix multiplication
    % h is input array, 
    % n is length of array to convolve 
    result = toeplitz([h, zeros(1, n-1)], [h(1), zeros(1, n-1)]);
end