function res = compute_joint_ci_ansatz(V,dim,x)
% V = isometry for joint channel
% dim = [B,E] with B = joint output dim. and E = joint environment dim.
% x = 4x1 variable vector

r = abs(x(1:3)); r = r/sum(r);
ep = abs(x(4));
if (ep>1)
    ep = 1/ep;
end
f0 = [1;0;0]; f1 = [0;1;0]; f2 = [0;0;1]; % computational basis for platypus input
e0 = [1;0]; e1 = [0;1]; % computational basis for qubit channel input
chi = sqrt(1-ep)*kron(f2,e0) + sqrt(ep)*kron(f1,e1);
rho = r(1)*(kron(f0,e0)*kron(f0,e0)') + r(2)*(kron(f0,e1)*kron(f0,e1)') + r(3)*(chi*chi');
sigma = V*rho*V';
res = real(VNent(TrX(sigma,2,dim))-VNent(TrX(sigma,1,dim)));
end