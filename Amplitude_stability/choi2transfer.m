function L = choi2transfer(T)
I = eye(4);
L = zeros(4);

for j=1:4
    L(:,j) = dm2bloch(applychan(T,bloch2dm(I(:,j)),'choi'));
end
end