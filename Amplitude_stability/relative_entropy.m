function f = pure_state_density(params)
    r = length(params);
    if (-1)^r == 1
        r = r/2;
    else
        disp('Not a valid even length real vector!')
    end
    % Extract real and imaginary parts
    realPart = params(1:r);
    imagPart = params(r + 1:2*r);
    % Construct the complex state vector
    psi = complex(realPart, imagPart);
    % Normalize the state vector
    f = psi / norm(psi);   
end

function S = quantum_relative_entropy(rho, sigma)
    % Check if rho and sigma are valid density matrices
    %if ~isequal(size(rho), size(sigma))
    %    error('Density matrices rho and sigma must have the same dimensions');
    %end
    %if ~ishermitian(rho) || ~ishermitian(sigma)
    %    error('Density matrices rho and sigma must be Hermitian');
    %end
    %if trace(rho) ~= 1 || trace(sigma) ~= 1
    %    error('Density matrices rho and sigma must have a trace of 1');
    %end
    
    % Calculate the logarithms of rho and sigma
    log_rho = logm(rho)/ log(2);
    log_sigma = logm(sigma)/ log(2);
    
    % Calculate the quantum relative entropy
    S = trace(rho * (log_rho - log_sigma));
end

