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
% Optimization 
opt = optimoptions('fminunc','disp','none');
opt_ps = optimoptions('particleswarm','UseParallel',true,'Display','none');

% Use a 50-50 erasure channel; set p to something else to test other erasure channels
lambda = 0.5;
E = erasure(2,lambda,'isom');
% compute coherent information of erasure channel
ci_E = max([1-2*lambda,0]);

% Use a depolarizing channel with p = 0.1893 (around hashing point)
% set p to some other value for other depolarizing channels
p = 0.1893; 
D = depol_chan(p,'isom');
% compute coherent information of depolarizing channel if not already known

% set channel dimensions, where b1, e1 corresponds to platypus and b2, e2
% to assisting channel
a1 = 2; b1 = 3; e1 = 3;
a2 = 2; b2 = 2; e2 = 4;
dim = [b1*b2,e1*e2];


% now look for amplification of coherent information
disp('Computing coherent information of joint channel...')

ED = kron(E,D);
    % need to reorder the systems from B1|E1|B2|E2 to B1|B2|E1|E2 because 
    % 'compute_priv_information.m' needs the output and environment systems
    % blocked together.
for k=1:size(ED,2)
        ED(:,k) = syspermute(ED(:,k),[1,3,2,4],[b1,e1,b2,e2]);
end
M = 10;
res = zeros(1,M);
res_x = cell(1, M);
for m=1:M
            obj = @(x) -compute_arbitrary_ci(ED,TrX(pure_state_density(x)'*pure_state_density(x), 2, [4,4]));
            x0 = rand(1,2*a1^2*a2^2);
            [x,f] = fminunc(obj,x0,opt);
            res(m) = -f;    
            res_x{m} = TrX(pure_state_density(x)'*pure_state_density(x), 2, [4,4]);
            disp(res(m))
            disp(res_x{m})
end
[mx, index]= max(res);
disp([mx,index])

function res = compute_arbitrary_ci(V,rho)
sigma = V*rho*V';
res = real(VNent(TrX(sigma,2,[6,12]))-VNent(TrX(sigma,1,[6,12])));
end

    

disp(['coherent information of depolarizing tensor erasure channel= ',num2str(mx)])






