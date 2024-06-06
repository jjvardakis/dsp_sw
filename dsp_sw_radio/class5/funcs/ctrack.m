function H = ctrack(A, T, I, P)
    % Function to calculate the open-loop transfer function
    
    % NCO transfer function (phase accumulator) (numerator z cancels with parasitic delay)
    num_nco = 2 * pi;
    den_nco = [1, -1];
    H_nco = tf(num_nco, den_nco, T);
    
    % PI Loop filter
    num_lf = [P, I * T - P];
    den_lf = [1, -1];
    H_lf = tf(num_lf, den_lf, T);
    
    % Open-loop transfer function
    H = A^2 * H_nco * H_lf;
end
