function binaryAngle = binAngle(angle, N)
    % Returns binary representation of angle between -pi to pi radians.
    % Each bit represents pi/2^N.
    binaryAngle = round((2^N-1)*(angle+pi/2)/pi);
    binaryAngle = dec2bin(binaryAngle, N);
end