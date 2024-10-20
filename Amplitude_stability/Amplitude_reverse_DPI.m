% generate 2 density matrices with r by r dimension; input is 4r^2
% dimensional real vector
function f = generate_two_density(params, mode)
    r = length(params)/4;
    % Extract real and imaginary parts
    realPart1 = params(1:r);
    imagPart1 = params(r + 1:2*r);
    realPart2 = params(2*r+1:3*r);
    imagPart2 = params(3*r+1:4*r);
    % Construct the complex state vector
    psi1 = complex(realPart1, imagPart1);
    psi2 = complex(realPart2, imagPart2);
    % Normalize the state vector
    psi1 = psi1 / norm(psi1);   
    psi2 = psi2 / norm(psi2);  
    rho1 = TrX(psi1'*psi1,2,[sqrt(r),sqrt(r)]);
    rho2 = TrX(psi2'*psi2,2,[sqrt(r),sqrt(r)]);
    switch mode
        case 1
            f = rho1;
        case 2
            f = rho2;
    end
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
% compute the advantange of amplitude damping with d-dimensional ancillary
% system given input rho; rho is 2d by 2d density matrix.
function f = relative_entropy_amplitude(gamma,rho,sigma)
    V = amp_damp(gamma,'isom');
    rho_BE = V*rho*V';
    sigma_BE = V*sigma*V';
    rho_B = TrX(rho_BE,2,[2,2]);
    sigma_B = TrX(sigma_BE,2,[2,2]);
    f = quantum_relative_entropy(rho_B, sigma_B);
end

%% 
% dimension of input is d=2
d = 2;
M = 50; % # number of optimizations 
res = zeros(1,M); % store the maximation for each random round
% Set optimization options
opt = optimoptions('fminunc','disp','none');
opt_ps = optimoptions('particleswarm','UseParallel',true,'Display','none');

% dimension of input is d=2
d = 2;
M = 50; % # number of optimizations 
res = zeros(1,M); % store the maximation for each random round
resolution = 101;
% Define the grid for gamma_1 and gamma_2
gamma_1 = linspace(0, 1, resolution);
gamma_2 = linspace(0, 1, resolution);
[Gamma1, Gamma2] = meshgrid(gamma_1, gamma_2);

disp(gamma_1)
%% 

% Preallocate the matrix to store the expansion values
expansion = zeros(size(Gamma1));

% Loop through the grid and find the maximum for each pair (gamma_1,gamma_2)
for i = 1:length(gamma_1) 
    for j = 1:length(gamma_2)
    % Apply the additional restriction gamma_1 > gamma_2
        if (Gamma1(i,j)> Gamma2(i,j)) && (Gamma2(i,j) > 0) && (Gamma1(i,j) < 1)
            disp(['Computing the inverse contraction coefficient of amplitute damping channel', ' with parameter=', num2str(Gamma1(i,j)), ',', num2str(Gamma2(i,j))])
            res = zeros(1,M);
            for m=1:M
               obj = @(x) relative_entropy_amplitude(Gamma1(i,j),generate_two_density(x,1),generate_two_density(x,2))/relative_entropy_amplitude(Gamma2(i,j),generate_two_density(x,1),generate_two_density(x,2));
               x0 = rand(1,4*d*d);
               % Perform the optimization
               [optimalParams, fval] = fminunc(obj,x0,opt);
               res(m)= fval;
            end
            expansion(i,j) = min(res);
            disp(['inverse contraction coefficient of amplitude damping channels= ', num2str(expansion(i,j))])
        elseif (Gamma1(i,j)<= Gamma2(i,j))
                 expansion(i,j) = 1;
        else
                 expansion(i,j) = 0;
        end
    end
end

