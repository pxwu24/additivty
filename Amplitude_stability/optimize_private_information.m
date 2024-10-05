function [pi,p,R] = optimize_private_information(V,dim,k,iter,mode,verb)

a = dim(1); b = dim(2); e = dim(3); r = dim(4);

disp_mode = 0;

if (nargin == 3)
    iter = 5;
    mode = 'fmin';
    verb = 'none';
elseif (nargin == 4)
    mode = 'fmin';
    verb = 'none';
    disp_mode = 1;
elseif (nargin == 5)
    verb = 'none';
    disp_mode = 1;
end

nx = 2*k*a*r;

switch mode
    case 'fmin'
        opt = optimoptions('fminunc','display',verb,'UseParallel',true);
    case 'psw'
        opt = optimoptions('particleswarm','display',verb,'UseParallel',true);
end



obj = @(x) -compute_priv_inf(x,V,dim,k);

res = zeros(1,iter);
res_x = zeros(nx,iter);
for l=1:iter
    if (disp_mode)
        disp(['Iteration: ',num2str(l)])
    end
    x0 = 100*randn(nx,1);
    switch mode
        case 'fmin'
            [x,f] = fminunc(obj,x0,opt);
        case 'psw'
            [x,f] = particleswarm(obj,nx,[],[],opt);
    end
    res(l) = -f;
    res_x(:,l) = x;
end
[mx,ix] = max(res);
pi = mx;
x = res_x(:,ix);

[p,R] = compute_ensemble(x,dim,k);
end