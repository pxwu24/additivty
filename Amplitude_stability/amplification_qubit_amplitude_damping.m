J = 50; % sample points in s-interval
s = linspace(0,1/2,J);
M = 10; % # optimizations for given value of s

% purifying dimension in the quantum state ensemble on which we evaluate 
% the private information in the full optimization.
% r = 1 means we're computing the coherent information
r = 1; 

opt = optimoptions('fminunc','disp','none');
opt_ps = optimoptions('particleswarm','UseParallel',true,'Display','none');

if (~exist('ci_vs'))
    disp('Computing coherent information of platypus channel...')
    vs_coherent_information
    fprintf('\n')
end

% Use a 50-50 amplitude damping channel; set p to something else to test 
% other amplitude damping channels
p = 0.4;
A = amp_damp(p,'isom');
% compute coherent information of amplitude damping channel
if (p>=0.5)
    ci_A = 0;
else
    ci_A = ci_amp_damp(p);
    disp(['ci(A) = ',num2str(ci_A)])
    fprintf('\n')
end

% set channel dimensions, where b1, e1 corresponds to platypus and b2, e2
% to assisting channel
a1 = 3; b1 = 3; e1 = 2;
a2 = 2; b2 = 2; e2 = 2;
dim = [b1*b2,e1*e2];

% now look for amplification of coherent information
disp('Computing coherent information of joint channel...')
for j=1:50
    disp(['s = ',num2str(s(j))])
    V = vs_channel(s(j),'isom');
    res = zeros(1,M);
    res_x = zeros(1,M);

    VA = kron(V,A);
    % need to reorder the systems from B1|E1|B2|E2 to B1|B2|E1|E2 because 
    % 'compute_priv_information.m' needs the output and environment systems
    % blocked together.
    for k=1:size(VA,2)
        VA(:,k) = syspermute(VA(:,k),[1,3,2,4],[b1,e1,b2,e2]);
    end
    
    % full optimization
    %     dim_full = [a1*a2,b1*b2,e1*e2,r]
    %     [ci,p,R] = optimize_private_information(WV,dim_full,a1*a2); % third argument is number of states in ensemble
    %     ci_random_full(j) = ci;
    %     disp(['ci(W x V) = ',num2str(ci),', ci(W) + ci(V) = ',num2str(ci_vs(j)+ci_V)])

    % optimization using ansatz eq.(28) in arXiv:2202.08377
    res = zeros(1,M);
    res_x = zeros(4,M);
    for m=1:M
        obj = @(x) -compute_joint_ci_ansatz(VA,dim,x);
        x0 = rand(1,4);
        [x,f] = fminunc(obj,x0,opt);
        res(m) = -f;
        res_x(:,m) = x;
    end
    ci_random_ansatz(j) = max(res);
    ci_random_ansatz_x(:,j) = x;
    disp(['ci(V x A) = ',num2str(ci_random_ansatz(j)),', ci(V) + ci(A) = ',num2str(ci_vs(j)+ci_A), ', Delta = ',num2str(ci_random_ansatz(j)-ci_vs(j)-ci_A)])
end

function ci = ci_amp_damp(p)
res = zeros(1,10);
A = amp_damp(p,'isom');
opt = optimoptions('fminunc','Display','none');
for j=1:10
    obj = @(x) -comp_ci_amp_damp(A,x);
    x0 = rand;
    [~,f] = fminunc(obj,x0,opt);
    res(j) = -f;
end
ci = max(res);
end

function res = comp_ci_amp_damp(A,x)
x = abs(x);
if (x>1)
    x = 1/x;
end
rho = diag([1-x,x]);
sigma = A*rho*A';
res = real(VNent(TrX(sigma,2,[2,2]))-VNent(TrX(sigma,1,[2,2])));
end