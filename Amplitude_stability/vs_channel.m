function res = vs_channel(S,mode)
f0 = [1;0;0]; f1 = [0;1;0]; f2 = [0;0;1];
e0 = [1;0]; e1 = [0;1];

V = [sqrt(S)*kron(f0,e0)+sqrt(1-S)*kron(f1,e1),...
    kron(f2,e0),...
    kron(f2,e1)];


if (nargin == 1)
    mode = 'isom';
end

switch mode
    case 'isom'
        res = V;
    case 'choi'
        res = chanconv(V,'isom','choi',[3,3]);
    case 'kraus'
        res = chanconv(V,'isom','kraus',[3,3]);
    otherwise
        error('Unknown channel representation!')
end
end