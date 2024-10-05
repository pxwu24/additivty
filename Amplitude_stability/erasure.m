function res = erasure(d,p,mode)
if (nargin < 3)
    mode = 'isom';
end
I = eye(d);
I_aug = [I;zeros(1,d)];
K{1,1} = sqrt(1-p)*I_aug;

for i=1:d
    A = zeros(d+1,d);
    A(d+1,i) = sqrt(p);
    K{i+1,1} = A;
end
switch mode
    case 'kraus'
        res = K;
    case 'isom'
        res = chanconv(K,'kraus','isom');
    case 'choi'
        res = chanconv(K,'kraus','choi');
end
end

