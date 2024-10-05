function res = compute_priv_inf(x,V,dim,k,flag)
if (nargin == 4)
    flag = 1;
end
b = dim(2); e = dim(3);
[p,R] = compute_ensemble(x,dim,k,flag);
rB_avg = zeros(b,b);
rE_avg = zeros(e,e);
ent_avg = 0;
for j=1:k
    op = V*R(:,:,j)*V';
    opB = TrX(op,2,[b,e]);
    opE = TrX(op,1,[b,e]);
    rB_avg = rB_avg + p(j)*opB;
    rE_avg = rE_avg + p(j)*opE;
    ent_avg = ent_avg + p(j)*real(VNent(opB)-VNent(opE));
end
res = real(VNent(rB_avg)-VNent(rE_avg)) - ent_avg;
end