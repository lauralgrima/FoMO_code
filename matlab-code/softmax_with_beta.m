function probabilities = softmax_with_beta(values, beta)
    % Softmax function with inverse temperature parameter beta
    %
    % Inputs:
    %   values: A vector of numerical values
    %   beta: Inverse temperature parameter (scalar)
    %
    % Output:
    %   probabilities: A vector of probabilities (same size as values)
    
    % Apply beta to the input values
    scaled_values = beta * values;
    
    % Subtract the maximum value for numerical stability
    shifted_values = scaled_values - max(scaled_values);
    
    % Calculate exponentials
    exp_values = exp(shifted_values);
    
    % Calculate sum of exponentials
    sum_exp = sum(exp_values);
    
    % Calculate probabilities
    probabilities = exp_values / sum_exp;
end