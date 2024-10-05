function v = dm2bloch(R)
I = eye(2);
[x,y,z] = paulis;

w0 = trace(I*R)/2;
w1 = trace(x*R)/2;
w2 = trace(y*R)/2;
w3 = trace(z*R)/2;

v = [w0;w1;w2;w3];
end