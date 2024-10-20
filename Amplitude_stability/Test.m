% find the s,t such that Q^{(1)}(N_{s,t}) = 0.
N = length(s_temp);
s_zero = [];
t_zero = [];
f_zero = [];
%% Tesing 
for i = 1:N
    if (f_temp(i)<0.000001) 
        s_zero = [s_zero; s_temp(i)];
        t_zero = [t_zero; t_temp(i)];
        f_zero = [f_zero; 0];
    end
end

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

M = 50; % # number of optimizations 
% Set optimization options
opt = optimoptions('fminunc','disp','none');
opt_ps = optimoptions('particleswarm','UseParallel',true,'Display','none');


% set channel dimensions, where b1, e1 corresponds to generalized platypus and b2, e2
% to assisting channel
a1 = 3; b1 = 3; e1 = 3;
a2 = 3; b2 = 3; e2 = 3;
dim_single = [b1,e1];
dim_two_copy = [b1*b2,e1*e2];


% Calculate S(B^2) - S(E^2) for a given non-diagonal rho(two copy case)
function res = compute_two_copies_vs_ci(V,rho)
sigma = V*rho*V';
res = real(VNent(TrX(sigma,2,[9,9]))-VNent(TrX(sigma,1,[9,9])));
end

% Calculate S(BB') - S(EE') for a given non-diagonal rho(tensor with erasure)
function res = compute_tensor_vs_ci(V,rho)
sigma = V*rho*V';
res = real(VNent(TrX(sigma,2,[9,9]))-VNent(TrX(sigma,1,[9,9])));
end

% Calculate S(B) - S(E) for a given non-diagonal rho
function res = compute_vs_ci(V,rho)
sigma = V*rho*V';
res = real(VNent(TrX(sigma,2,[3,3]))-VNent(TrX(sigma,1,[3,3])));
end

% specify s,t and other parameters
d= 3;
p = 0.5;
disp(length(s_zero))

for u = 1: length(s_zero)
         % specify the channel isometry
         V = Generalized_vs_channel(s_zero(u),t_zero(u),'isom');
         E = erasure(2,p,'isom');
         VE = kron(V,E);
         % need to reorder the systems from B1|E1|B2|E2 to B1|B2|E1|E2 because 
    % 'compute_priv_information.m' needs the output and environment systems
    % blocked together.
         for k=1:size(VE,2)
             VE(:,k) = syspermute(VE(:,k),[1,3,2,4],[b1,e1,b2,e2]);
         end
% V_2: 81 * 9 isometry 
         V_2 = kron(V,V);
% change B_1E_1B_2E_2 to B_1B_2E_1E_2
         for k=1:size(V_2,2)
             V_2(:,k) = syspermute(V_2(:,k),[1,3,2,4],[b1,e1,b2,e2]);
         end


% specify the value and cell to store the coherent information
% res_general is for non-diagonal optimization for single copy
         res_general = zeros(1,M);
         res_general_rho = cell(1, M);
% res_general_two_copy is for non-diagonal optimization for two copy
         res_general_two_copy = zeros(1,M);
         res_general_rho_two_copy = cell(1, M);

% res_general_erasure is for non-diagonal optimization for tensor product
% with erasure
         res_general_erasure = zeros(1,M);
         res_general_rho_erasure = cell(1, M);


         disp('Computing coherent information of generalized Platypus channel...')
         for m=1:M
            obj = @(x) -compute_vs_ci(V,TrX(pure_state_density(x)'*pure_state_density(x), 2, [3,3]));
            x0 = rand(1,2*9);
            [x,f] = fminunc(obj,x0,opt);
            res_general(m) = -f;    
            res_general_rho{m} = TrX(pure_state_density(x)'*pure_state_density(x), 2, [3,3]);
         end      

         disp('Computing coherent information of two copies of generalized Platypus channel...')
         for m=1:M
            obj = @(x) -compute_two_copies_vs_ci(V_2,TrX(pure_state_density(x)'*pure_state_density(x), 2, [9,9]));
            x0 = rand(1,2*81);
         % Perform the optimization
            [x, f] = fminunc(obj,x0,opt);
         % Desired density:
            res_general_two_copy(m) = -f;
            res_general_rho_two_copy{m} = TrX(pure_state_density(x)'*pure_state_density(x), 2, [9,9]);
         end

         disp('Computing coherent information of generalized Platypus channel tensor with erasure channel...')
         for m=1:M
            obj = @(x) -compute_tensor_vs_ci(VE,TrX(pure_state_density(x)'*pure_state_density(x), 2, [6,6]));
            x0 = rand(1,2*36);
            [x,f] = fminunc(obj,x0,opt);
            res_general_erasure(m) = -f;
            res_general_rho_erasure{m} = TrX(pure_state_density(x)'*pure_state_density(x), 2, [6,6]);
         end
         [ci_single, index_single]= max(res_general);
         [ci_two_copy, index_two_copy]= max(res_general_two_copy);
         [ci_erasure, index_erasure]= max(res_general_erasure);

         disp(['The parameters of generalized Platypus channel are: ', 's= ', num2str(s_zero(u)),', t= ', num2str(t_zero(u))])
         disp(['coherent information of generalized Platypus channel = ', num2str(ci_single)])
         disp(['coherent information of two-copy generalized Platypus channel = ', num2str(ci_two_copy)])
         disp(['coherent information of generalized Platypus channel tensor with erasure= ', num2str(ci_erasure)])
end
%disp(res_general_rho{index_gen})

