function [p,R,Ravg] = compute_ensemble(x,dim,k,flag)
if (nargin == 3)
    flag = 1;
end
if (flag)
    psi(:,1) = x(1:length(x)/2) + 1i*x(length(x)/2+1:end);
else
    psi(:,1) = x;
end
psi = psi/norm(psi);

a = dim(1); r = dim(4);

R = zeros(a,a,k);
p = zeros(1,k);
Ravg = zeros(a,a);
I = eye(k);
for j=1:k
    phi = kron(I(:,j)',eye(a*r))*psi;
    if (r>1)
        rho = TrX(phi,2,[a,r]);
    else
        rho = phi*phi';
    end
    p(j) = trace(rho);
    R(:,:,j) = rho/trace(rho);
    Ravg = Ravg + p(j)*R(:,:,j);
end
end