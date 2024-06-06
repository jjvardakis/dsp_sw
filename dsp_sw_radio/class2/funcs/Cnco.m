function nco_gen = Cnco(acc_size, lut_in_size, lut_out_size, dither, state)
    % state = accumulator output, lut input
    
    % Initialize the state if not provided
    if nargin < 5
        state = [0, 0];
    end
    
    % Function handle for the NCO generator
    nco_gen = @nco_generator;

    % Define the generator function
    function result = nco_generator(fcw, pcw, dither)
        angle = state(2) / 2^lut_in_size * 2 * pi;
        magnitude = 2^(lut_out_size-1);
        result = floor(magnitude * cos(angle)) + 1i * floor(magnitude * sin(angle));

        if nargin > 0
            if nargin < 2
                pcw = 0;
            end

            if ~isempty(dither)
                dith = round(2^(dither-1) * randn(1));
            else
                dith = 0;
            end

            % Update phase accumulator
            state(1) = mod(fcw + state(1) + dith, 2^acc_size);

            % Add phase control word and truncate phase
            state(2) = floor(mod((pcw + state(1)), 2^acc_size) / 2^(acc_size-lut_in_size));
        end
    end
end

