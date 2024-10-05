% generate r dimensional pure state from 2r real vectors.
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

% compute the advantange of amplitude damping with d-dimensional ancillary
% system given input rho
function f = compute_advantage(gamma,d,rho)
    V = amp_damp(gamma,'isom');
    V_amp = kron(eye(d), V);
    sigma_VBE = V_amp*rho*V_amp';
    sigma_VB = TrX(sigma_VBE,2,[2*d,2]);
    sigma_VE = TrX(sigma_VBE,2,[d,2,2]);
    sigma_BE = TrX(sigma_VBE,1,[d,4]); 
    % f1 is S(B)- S(E); f2 is S(BV) - S(EV)
    f1 = real(VNent(TrX(sigma_BE,2,[2,2]))-VNent(TrX(sigma_BE,1,[2,2])));
    f2 = real(VNent(sigma_VB)-VNent(sigma_VE));
    f = f1 - f2;
end
% dimension of ancillary system
d = 4;
psi = pure_state_density(rand(1, 2*(2*d)^2));
rho_AV = TrX(psi'*psi, 2, [2*d,2*d]);

M = 10; % # number of optimizations 
% Set optimization options
opt = optimoptions('fminunc','disp','none');
opt_ps = optimoptions('particleswarm','UseParallel',true,'Display','none');

disp('Computing the minimum ratio of advantange of two amplitude damping channels')
    res = zeros(1,M);
    for m=1:M
        obj = @(x) compute_advantage(0.21,d,TrX(pure_state_density(x)'*pure_state_density(x), 2, [2*d,2*d]))/compute_advantage(0.2,d,TrX(pure_state_density(x)'*pure_state_density(x), 2, [2*d,2*d]));
        x0 = rand(1,2*(2*d)^2);
        % Perform the optimization
        [optimalParams, fval] = fminunc(obj,x0,opt);

        % Extract optimized real and imaginary parts
        optimalRealPart = optimalParams(1:(2*d)^2);
        optimalImagPart = optimalParams((2*d)^2+1:2*(2*d)^2);

        % Construct the optimized complex state vector
        optimalPsi = complex(optimalRealPart, optimalImagPart);

        % Normalize the optimized state vector
        optimalPsi = optimalPsi / norm(optimalPsi);

        % Desired density:
        %res_state(m)= TrX(optimalPsi'*optimalPsi,2,[d,d])
        res(m) = fval;
    end


disp(['minimum ratio of advantange of two amplitude damping channels= ', num2str(min(res))])