function Ravg = compute_average_state(p,R)
assert(length(p)==size(R,3),'Lists of probabilities and states are not of same length.')

Ravg = zeros(size(R,[1,2]));
for j=1:length(p)
    Ravg = Ravg + p(j)*R(:,:,j);
end
end