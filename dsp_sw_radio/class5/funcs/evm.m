function result = evm(x, ref)
    % computes error vector magnitude for x and ref 
    % where x is the samples at the decision locations and
    % ref is the reference correct decision
    result = std(x - ref) / std(ref);
end