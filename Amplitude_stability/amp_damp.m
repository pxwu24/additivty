function res = amp_damp(p,mode)
% Computes the Choi state of the qubit amplitude damping channel with
% parameter p.

K1 = [1,0; 0, sqrt(1-p)];
K2 = [0, sqrt(p); 0, 0];

I = eye(2);
phi = [1;0;0;1];
Phi = phi*phi';

if (nargin<2)
    mode = 'choi';
end

T = tensor(I,K1)*Phi*tensor(I,K1)' + tensor(I,K2)*Phi*tensor(I,K2)';
K = {K1,K2};

switch mode
    case 'choi'
        res = T;
    case 'kraus'
        res = K;
    case 'isom'
        res = chanconv(K,'kraus','isom',[2,2]);
end
end