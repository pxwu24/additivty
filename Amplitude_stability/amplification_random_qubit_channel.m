J = 50; % sample points in s-interval
s = linspace(0,1/2,J);
M = 20; % # optimizations for given value of s

% purifying dimension in the quantum state ensemble on which we evaluate 
% the private information.
% r = 1 means we're computing the coherent information
r = 1; 

opt = optimoptions('fminunc','disp','none');
opt_ps = optimoptions('particleswarm','UseParallel',true,'Display','none');

if (~exist('ci_vs'))
    disp('Computing coherent information of platypus channel...')
    vs_coherent_information
    fprintf('\n')
end

flag = false;
while (~flag)
    l = -1+2*rand(1,3);
    t = -1+2*rand(3,1);
    L = [1,zeros(1,3);t,diag(l)];
    flag = IsPSD(transfer2choi(L));
end

T = transfer2choi(L);

assert(IsPSD(T),'Channel is not CP')
assert(norm(TrX(T,2,[2,2])-eye(2))<1e-10,'Channel is not TP')

% First determine the "hashing point" of the channel family
% (1-x)*Id + x*T
p = linspace(0,1,50);
Phi = [1,0,0,1;0,0,0,0;0,0,0,0;1,0,0,1];
ci_a = zeros(1,50);
disp('Determining hashing point of assisting channel...')
for j=1:50
    D = (1-p(j))*Phi + p(j)*T;
    A = chanconv(D,'choi','isom');
    e2 = size(A,1)/2;
    ci_a(j) = optimize_private_information(A,[2,2,e2,1],2);
    disp(['j = ',num2str(j),', p = ',num2str(p(j)),', ci = ',num2str(ci_a(j))])
    if (ci_a(j) < 0)
        ci_a(j) = 0;
        j0 = j;
        break;
    end
end
fprintf('\n')

% set channel dimensions, where b1, e1 corresponds to platypus and b2, e2
% to assisting channel
a1 = 3; b1 = 3; e1 = 2;
a2 = 2; b2 = 2;
dim = [b1*b2,e1*e2];

%%
% now look for amplification of coherent information
disp('Computing coherent information of joint channel...')
for j=50
    disp(['s = ',num2str(s(j))])
    V = vs_channel(s(j),'isom');
    res = zeros(1,M);
    res_x = zeros(1,M);

    VA = kron(V,A);
    % reordering the systems from B1|E1|B2|E2 to B1|B2|E1|E2 so that 
    % output and environment systems are blocked together.
    for k=1:size(VA,2)
        VA(:,k) = syspermute(VA(:,k),[1,3,2,4],[b1,e1,b2,e2]);
    end
    
    % full optimization
    %     dim_full = [a1*a2,b1*b2,e1*e2,r];
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
    disp(['ci(V x A) = ',num2str(ci_random_ansatz(j)),', ci(V) + ci(A) = ',num2str(ci_vs(j)+ci_a(j0)), ', Delta = ',num2str(ci_random_ansatz(j)-ci_vs(j)-ci_a(j0))])
end