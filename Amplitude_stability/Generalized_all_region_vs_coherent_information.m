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


% Define the resolution of the grid
resolution = 100; % Increase for finer grids/ number of refinement
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



% res_general is for non-diagonal optimization
res_general = zeros(1,M);
res_general_rho = cell(1, M);

% Create arrays to store valid (s, t) pairs and corresponding Q^(1)(N_{s,t})
s_values = [];
t_values = [];
f_values = [];

% Loop over possible values of s and t
for i = 1:resolution
    for j = 1:resolution
        s = i / resolution;  % Generate values of s between 0 and 1
        t = j / resolution;  % Generate values of t between 0 and 1
        
        % Apply the conditions s + t <= 1 and t < 1/2
        if (s + t <= 1) && (t < 0.5)
            s_values = [s_values; s];   % Store valid s
            t_values = [t_values; t];   % Store valid t
            V = Generalized_vs_channel(s,t,'isom'); % isomotry of generalized Platypus
            for m=1:M
            obj = @(x) -compute_vs_ci(V,TrX(pure_state_density(x)'*pure_state_density(x), 2, [3,3]));
            x0 = rand(1,2*3*3);
            [x,f] = fminunc(obj,x0,opt);
            res_general(m) = -f;    
            res_general_rho{m} = TrX(pure_state_density(x)'*pure_state_density(x), 2, [3,3]);
            end      
            [mx_gen, index_gen]= max(res_general);
            disp(['The parameters of generalized Platypus channel are: ', 's= ', num2str(s),', t= ', num2str(t)])
            disp(['coherent information of generalized Platypus channel = ', num2str(mx_gen)])
            disp(res_general_rho{index_gen})
            f_values = [f_values; mx_gen]; % Store the function value Q^(1)(N_{s,t})
        end
    end
end

% Create a scatter plot of the points in the region
scatter3(s_values, t_values, f_values, 36, f_values, 'filled');

% Label axes
xlabel('s');
ylabel('t');
zlabel('Q^(1)(N_{s,t})');
title('Plot of Q^(1)(N_{s,t}) in the region s+t \leq 1 and t < 1/2');

% Add a color bar to indicate the function value
colorbar;

% Set view for better visualization
view(3);

% Calculate S(B) - S(E) for a given non-diagonal rho
function res = compute_vs_ci(V,rho)
sigma = V*rho*V';
res = real(VNent(TrX(sigma,2,[3,3]))-VNent(TrX(sigma,1,[3,3])));
end

