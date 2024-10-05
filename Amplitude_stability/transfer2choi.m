function T = transfer2choi(L)
e0 = [1;0]; e1 = [0;1];
E00 = e0*e0'; E01 = e0*e1';
E10 = e1*e0'; E11 = e1*e1';

E = {E00,E01,E10,E11};
T = zeros(4);

for j=1:4
    w = L*dm2bloch(E{j});
    T = T + kron(E{j},bloch2dm(w));
end
end