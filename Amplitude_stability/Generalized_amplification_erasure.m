L = 30; % number of partition of the region of (s,t)
M = 5; % # optimizations for given value of s
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
p = 0.5;
E = erasure(2,p,'isom');
% compute coherent information of erasure channel
ci_E = max([1-2*p,0]);

% set channel dimensions, where b1, e1 corresponds to generalized platypus and b2, e2
% to assisting channel
a1 = 3; b1 = 3; e1 = 3;
a2 = 2; b2 = 3; e2 = 3;
dim = [b1*b2,e1*e2];



% now look for amplification of coherent information
disp('Computing coherent information of joint channel...')


V = Generalized_vs_channel(0.25, 0.25, 'isom');
VE = kron(V,E);
    % need to reorder the systems from B1|E1|B2|E2 to B1|B2|E1|E2 because 
    % 'compute_priv_information.m' needs the output and environment systems
    % blocked together.
for k=1:size(VE,2)
        VE(:,k) = syspermute(VE(:,k),[1,3,2,4],[b1,e1,b2,e2]);
end
res = zeros(1,M);
res_x = zeros(1,M);
for m=1:M
            obj = @(x) -compute_arbitrary_ci(VE,TrX(pure_state_density(x)'*pure_state_density(x), 2, [6,6]));
            x0 = rand(1,2*a1^2*a2^2);
            [x,f] = fminunc(obj,x0,opt);
            res(m) = -f;    
end
mx= max(res);

function res = compute_arbitrary_ci(V,rho)
sigma = V*rho*V';
res = real(VNent(TrX(sigma,2,[9,9]))-VNent(TrX(sigma,1,[9,9])));
end

    

disp(['coherent information of generalized Platypus channel tensor 0.5 erasure channel= ',num2str(mx)])


