figure;
circle = linspace(0, 2*pi, 100);
plot(sin(circle), cos(circle), 'r--');
hold on
N = 8;
step = pi / N;
offset = step / 2;
for n = 0:N-1
    quiver(0, 0, sin(n*step + offset), cos(n*step + offset), 'b', 'LineWidth', 1.5, 'ShowArrowHead','on');
    hold on
end
axis equal;
xlabel('I');
ylabel('Q');
title('CORDIC Rotation Angles - Ideal');


% Simple 3 iteration CORDIC demonstration
iterations = 3;
I = 1;
Q = 0;

fprintf('%-10s %-12s %-12s %-10s\n', 'Sequence', 'I', 'Q', 'Angle (deg)');
for binWord = 0:2^iterations-1
    sequence = dec2bin(binWord, iterations);
    [Ir, Qr] = cordicSequence(I, Q, sequence);
    angle_deg = atan2d(Qr, Ir);
    fprintf('%-10s %-12.2f %-12.2f %-10.2f\n', sequence, Ir, Qr, angle_deg);
end

%% Count through a complete circle to compare with ideal plot
N = 4; % number of iterations +1 and number of bits in sequence

figure;

for count = 0:2^N-1
    sequence = dec2bin(count, N);
    
    % plotting just the rotation vector, so start from 1 + j0, or -1 + j0 angles beyond +/-90°
    if (sequence(1) == '1')
        I1 = -1; % rotate vector 180 degrees for full 360° rotation
    else
        I1 = 1;
    end
    Q1 = 0;
        
    [Ir, Qr] = cordicSequence(I1, Q1, sequence(2:end));
    
    % plot result
    quiver(0, 0, Ir, Qr, 'Color', 'b', 'MaxHeadSize', 0.5, 'LineWidth', 1.5);
    hold on;
end

circle = linspace(0, 2*pi, 100);
cordicGain = 1.647;
plot(cordicGain * sin(circle), cordicGain * cos(circle), 'r--');    
axis equal;
xlabel('I');
ylabel('Q');
title('Cordic Rotation Angles - Actual');
hold off;


%% 5 iteration CORDIC
iterations = 5;
result = zeros(1, 2^iterations); % preallocate array to store result
I = 1;
Q = 0;

fprintf('%-10s %-12s %-12s %-10s\n', 'Sequence', 'I', 'Q', 'Angle (deg)');
for binWord = 0:2^iterations-1
    sequence = dec2bin(binWord, iterations);
    [Ir, Qr] = cordicSequence(I, Q, sequence);
    result(binWord+1) = complex(Ir, Qr);
    angle_deg = atan2d(Qr, Ir);
    fprintf('%-10s %-12.2f %-12.2f %-10.2f\n', sequence, Ir, Qr, angle_deg);
end