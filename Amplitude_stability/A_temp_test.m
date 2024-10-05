L = 30; % number of partition of the region of (s,t)
M = 50; % # optimizations for given value of s


opt = optimoptions('fminunc','disp','none');
% test s=t= 0.25
V = Generalized_vs_channel(0.25, 0.25, 'isom');
res = zeros(1,M);
        for m=1:M
            obj = @(x) -compute_vs_ci(V,x);
            x0 = rand(1,2);
            [x,f] = fminunc(obj,x0,opt);
            res(m) = -f;    
        end
        mx = max(res);
disp(['coherent information of a single generalized Platypus channel=',num2str(mx)])



function res = compute_vs_ci(V,x)
x = abs(x);
if (sum(x)>1)
    x = x/sum(x);
end
rho = diag([x(1),x(2),1-x(1)-x(2)]);
sigma = V*rho*V';
res = real(VNent(TrX(sigma,2,[3,3]))-VNent(TrX(sigma,1,[3,3])));
end