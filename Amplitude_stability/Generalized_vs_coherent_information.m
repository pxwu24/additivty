L = 30; % number of partition of the region of (s,t)
M = 50; % # optimizations for given value of s


opt = optimoptions('fminunc','disp','none');

V = Generalized_vs_channel(0.25, 0.25, 'isom');
res = zeros(1,M);
        res_x = zeros(1,M);
        for m=1:M
            obj = @(x) -compute_vs_ci(V,x);
            x0 = rand(1,2);
            [x,f] = fminunc(obj,x0,opt);
            res(m) = -f;    
        end 
        mx = max(res);


for j=1:L
    for k= 1:(L-j+1)
        V = Generalized_vs_channel((j-1)/L, k/L, 'isom');
        res = zeros(1,M);
        res_x = zeros(1,M);
        for m=1:M
            obj = @(x) -compute_vs_ci(V,x);
            x0 = rand(1,2);
            [x,f] = fminunc(obj,x0,opt);
            res(m) = -f;    
        end
        mx = max(res);
        disp(['s = ',num2str((j-1)/L),',t = ', num2str(k/L),',ci = ',num2str(mx)])
    end
end

function res = compute_vs_ci(V,x)
x = abs(x);
if (sum(x)>1)
    x = x/sum(x);
end
rho = diag([x(1),x(2),1-x(1)-x(2)]);
sigma = V*rho*V';
res = real(VNent(TrX(sigma,2,[3,3]))-VNent(TrX(sigma,1,[3,3])));
end