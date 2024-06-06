function result = qam(num_symb, order)
    % create random QAM data symbols
    qam_data = randi([0, order-1], 1, num_symb);

    % order: 4=QPSK, 16=16AM, etc
    result = qammod(qam_data, order);
end