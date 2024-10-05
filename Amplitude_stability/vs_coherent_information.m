J = 50; % sample points in s-interval
s = linspace(0,1/2,J);
M = 10; % # optimizations for given value of s

opt = optimoptions('fminunc','disp','none');

ci_vs = zeros(1,J);
ci_vs_x = zeros(1,J);

for j=1:J
    V = vs_channel(s(j),'isom');
    res = zeros(1,M);
    res_x = zeros(1,M);
    
    for m=1:M
        obj = @(x) -compute_vs_ci(V,x);
        x0 = rand;
        [x,f] = fminunc(obj,x0,opt);
        res(m) = -f;
        res_x(m) = x;
    end
    [mx,ix] = max(res);
    disp(['s = ',num2str(s(j)),', ci = ',num2str(mx)])
    ci_vs(j) = mx;
    ci_vs_x(j) = x;
end

function res = compute_vs_ci(V,x)
x = abs(x);
if (x>1)
    x = 1/x;
end
rho = diag([1-x,0,x]);
sigma = V*rho*V';
res = real(VNent(TrX(sigma,2,[3,2]))-VNent(TrX(sigma,1,[3,2])));
end