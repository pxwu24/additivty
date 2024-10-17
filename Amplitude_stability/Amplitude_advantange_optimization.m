% generate r dimensional pure state from 2r real vectors. params is a 2r
% real vector
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
% system given input rho; rho is 2d by 2d density matrix.
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
d = 2;
psi = pure_state_density(rand(1, 2*(2*d)^2));
rho_AV = TrX(psi'*psi, 2, [2*d,2*d]);

M = 20; % # number of optimizations 
% Set optimization options
opt = optimoptions('fminunc','disp','none');
opt_ps = optimoptions('particleswarm','UseParallel',true,'Display','none');

n = 50;
gamma = [];
eta = [];
ratio_array = [];
opt_rho = cell(1,M);
disp('Computing the minimum ratio of advantange of two amplitude damping channels')
for k=1:(n/2-1)
 for u = 0:(k-1)   
    res = zeros(1,M);
    for m=1:M
        obj = @(x) compute_advantage(k/n,d,TrX(pure_state_density(x)'*pure_state_density(x), 2, [2*d,2*d]))/compute_advantage(u/n,d,TrX(pure_state_density(x)'*pure_state_density(x), 2, [2*d,2*d]));
        x0 = rand(1,2*(2*d)^2);
        % Perform the optimization
        [optimalParams, fval] = fminunc(obj,x0,opt);
        % Desired density:
        opt_rho{m} = TrX(pure_state_density(optimalParams)'*pure_state_density(optimalParams), 2, [2*d,2*d]);
        % Optimal value:
        res(m) = fval;
    end

    [temp_ratio,Index] = min(res);
    gamma = [gamma,k/n];
    eta = [eta,u/n];
    ratio_array = [ratio_array, temp_ratio];
    rho_V = TrX(opt_rho{Index},2,[d,2]); 
    rho_A = TrX(opt_rho{Index},1,[d,2]);
    disp(['minimum ratio of informational advantange of two amplitude damping channels with parameter', ' ', num2str(k/n),' and ', num2str(u/n), ' is: ', num2str(min(res))])
    disp('Optimizer rho_VA minus rho_V tensor rho_A is:')
    disp(opt_rho{Index} - kron(rho_V, rho_A))
 end
end
% eta is \gamma_1 and gamma is \gamma_2 so that 0 < \gamma_1 < \gamma_2 < \frac{1}{2}
% Convert x, y, z to matrices suitable for surf
[Xq, Yq] = meshgrid(unique(eta), unique(gamma));
Zq = griddata(eta, gamma, ratio_array, Xq, Yq);

% Plotting
surf(Xq, Yq, Zq);
xlabel('$\gamma_1$','Interpreter','latex');
ylabel('$\gamma_2$','Interpreter','latex');
zlabel('$R(\gamma_1,\gamma_2)$','Interpreter','latex');
colorbar;
title('Surface Plot of $R(\gamma_1,\gamma_2)$ with $0 < \gamma_1 < \gamma_2 < \frac{1}{2}$','Interpreter','latex');