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

M = 10; % # number of optimizations 
% Set optimization options
opt = optimoptions('fminunc','disp','none');
opt_ps = optimoptions('particleswarm','UseParallel',true,'Display','none');


% set channel dimensions, where b1, e1 corresponds to generalized platypus and b2, e2
% to assisting channel
a1 = 3; b1 = 3; e1 = 3;
a2 = 3; b2 = 3; e2 = 3;
dim = [b1*b2,e1*e2];


% Calculate S(B) - S(E) for a given rho
function res = compute_two_copies_vs_ci(V,rho)
sigma = V*rho*V';
res = real(VNent(TrX(sigma,2,[9,9]))-VNent(TrX(sigma,1,[9,9])));
end

% Try simple s,t, d= 3
% two copies of generalized Platypus
d= 9;
s= 0.25;
t= 0.25;
V = Generalized_vs_channel(s,t,'isom');

% V_2: 81 * 9 isometry 
V_2 = kron(V,V);

% change B_1E_1B_2E_2 to B_1B_2E_1E_2
for k=1:size(V_2,2)
        V_2(:,k) = syspermute(V_2(:,k),[1,3,2,4],[b1,e1,b2,e2]);
end
psi = pure_state_density(rand(1, 2*d^2));
rho = TrX(psi'*psi, 2, dim);




% now look for coherent information of two copies
disp('Computing coherent information of two copies of the channel...')
    res = zeros(1,M);
    
    for m=1:M
        obj = @(x) -compute_two_copies_vs_ci(V_2,TrX(pure_state_density(x)'*pure_state_density(x), 2, dim));
        x0 = rand(1,2*d^2);
        % Perform the optimization
        [optimalParams, fval] = fminunc(obj,x0,opt);

        % Extract optimized real and imaginary parts
        optimalRealPart = optimalParams(1:d^2);
        optimalImagPart = optimalParams(d^2+1:2*d^2);

        % Construct the optimized complex state vector
        optimalPsi = complex(optimalRealPart, optimalImagPart);

        % Normalize the optimized state vector
        optimalPsi = optimalPsi / norm(optimalPsi);

        % Desired density:
        %res_state(m)= TrX(optimalPsi'*optimalPsi,2,[d,d])
        res(m) = fval;
    end


disp(['coherent information of two copies of generalized Platypus channel= ', num2str(max(-res))])

