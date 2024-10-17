% generate r dimensional pure state from 2r real vectors.
function f = pure_state_density(params)
    r = length(params);
    if (-1)^r == 1
        r = r/2;
    else
        disp('Not a valid even length real vector!')
    end
    % Extract real and imaginary parts
    realPart = params(1:r);
    imagPart = params(r + 1:2*r);
    % Construct the complex state vector
    psi = complex(realPart, imagPart);
    % Normalize the state vector
    f = psi / norm(psi);   
end



L = 30; % number of partition of the region of (s,t)
M = 100; % # optimizations for given value of s,t


opt = optimoptions('fminunc','disp','none');
% Choice of the parameters s and t.
% while true
%     s = rand(); % Generate a random number between 0 and 1
%    t = rand()/2; % Generate a random number between 0 and 1
%    if s + t <= 1
%        break; % If the condition is satisfied, exit the loop
%    end
% end
s = 0.6;
t = 0.3;
V = Generalized_vs_channel(s,t,'isom');

% res is for diagonal optimization

res = zeros(1,M);
res_rho = cell(1, M);
for m=1:M
       obj = @(x) -compute_vs_ci_diag(V,x);
       x0 = rand(1,2);
       [x,f] = fminunc(obj,x0,opt);
       res(m) = -f;
       res_rho{m} = diag([x(1),x(2),1-x(1)-x(2)]);
end 
[mx_diag, index_diag]= max(res);

% res_general is for non-diagonal optimization
res_general = zeros(1,M);
res_general_rho = cell(1, M);
for m=1:M
            obj = @(x) -compute_vs_ci(V,TrX(pure_state_density(x)'*pure_state_density(x), 2, [3,3]));
            x0 = rand(1,2*3*3);
            [x,f] = fminunc(obj,x0,opt);
            res_general(m) = -f;    
            res_general_rho{m} = TrX(pure_state_density(x)'*pure_state_density(x), 2, [3,3]);
end      
[mx_gen, index_gen]= max(res_general);


%disp(['coherent information of generalized Platypus channel restricted on diagonal operator= ', num2str(mx_diag)])
%disp(res_rho{index_diag})
disp(['The parameters of generalized Platypus channel are: ', 's= ', num2str(s),' t= ', num2str(t)])
disp(['coherent information of generalized Platypus channel = ', num2str(mx_gen)])
disp(res_general_rho{index_gen})
%for j=1:L
%    for k= 1:(L-j+1)
%        V = Generalized_vs_channel((j-1)/L, k/L, 'isom');
%        res = zeros(1,M);
%        res_x = zeros(1,M);
%        for m=1:M
%            obj = @(x) -compute_vs_ci(V,x);
%            x0 = rand(1,2);
%            [x,f] = fminunc(obj,x0,opt);
%            res(m) = -f;    
%        end
%        mx = max(res);
%        disp(['s = ',num2str((j-1)/L),',t = ', num2str(k/L),',ci = ',num2str(mx)])
%    end
% end


% Calculate S(B) - S(E) for a given non-diagonal rho
function res = compute_vs_ci(V,rho)
sigma = V*rho*V';
res = real(VNent(TrX(sigma,2,[3,3]))-VNent(TrX(sigma,1,[3,3])));
end

% Calculate S(B) - S(E) for a given diagonal(x)
function res = compute_vs_ci_diag(V,x)
x = abs(x);
if (sum(x)>1)
    x = x/sum(x);
end
rho = diag([x(1),x(2),1-x(1)-x(2)]);
sigma = V*rho*V';
res = real(VNent(TrX(sigma,2,[3,3]))-VNent(TrX(sigma,1,[3,3])));
end