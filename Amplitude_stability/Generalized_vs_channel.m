function res = Generalized_vs_channel(S,T, mode)
f0 = [1;0;0]; f1 = [0;1;0]; f2 = [0;0;1];
e0 = [1;0;0]; e1 = [0;1;0]; e2 = [0;0;1];

V = [sqrt(S)*kron(f0,e0)+sqrt(1-S-T)*kron(f1,e1)+sqrt(T)*kron(f2,e2),...
    kron(f2,e0),...
    kron(f2,e1)];


if (nargin == 1)
    mode = 'isom';
end

switch mode
    case 'isom'
        res = V;
    otherwise
        error('Unknown channel representation!')
end
end