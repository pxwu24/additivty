function R = bloch2dm(v)
I = eye(2);
[x,y,z] = paulis;

R = v(1)*I+v(2)*x+v(3)*y+v(4)*z;
end