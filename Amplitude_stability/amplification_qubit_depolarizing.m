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

% Use a depolarizing channel with p = 0.1893 (around hashing point)
% set p to some other value for other depolarizing channels
p = 0.1893; 
D = depol_chan(p,'isom');
% compute coherent information of depolarizing channel if not already known
if (p>=0.25)
    ci_D = 0;
else
    ci_D = max([0,ci_depol(p)]);
    disp(['ci(D) = ',num2str(ci_D)])
    fprintf('\n')
end

% set channel dimensions, where b1, e1 corresponds to platypus and b2, e2
% to assisting channel
a1 = 3; b1 = 3; e1 = 2;
a2 = 2; b2 = 2; 
if (p>0)
    e2 = 4;
else
    e2 = 1;
end
dim = [b1*b2,e1*e2];

% now look for amplification of coherent information
disp('Computing coherent information of joint channel...')
for j=1:50
    disp(['s = ',num2str(s(j))])
    V = vs_channel(s(j),'isom');
    res = zeros(1,M);
    res_x = zeros(1,M);

    VD = kron(V,D);
    % need to reorder the systems from B1|E1|B2|E2 to B1|B2|E1|E2 because 
    % 'compute_priv_information.m' needs the output and environment systems
    % blocked together.
    for k=1:size(VD,2)
        VD(:,k) = syspermute(VD(:,k),[1,3,2,4],[b1,e1,b2,e2]);
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
        obj = @(x) -compute_joint_ci_ansatz(VD,dim,x);
        x0 = rand(1,4);
        [x,f] = fminunc(obj,x0,opt);
        res(m) = -f;
        res_x(:,m) = x;
    end
    ci_random_ansatz(j) = max(res);
    ci_random_ansatz_x(:,j) = x;
    disp(['ci(V x D) = ',num2str(ci_random_ansatz(j)),', ci(V) + ci(D) = ',num2str(ci_vs(j)+ci_D), ', Delta = ',num2str(ci_random_ansatz(j)-ci_vs(j)-ci_D)])
end

function ci = ci_depol(p)
if (p==0)
    ci = 1;
else
    rho = eye(2)/2;
    D = depol_chan(p,'isom');
    sigma = D*rho*D';
    ci = real(VNent(TrX(sigma,2,[2,4]))-VNent(TrX(sigma,1,[2,4])));
end
end