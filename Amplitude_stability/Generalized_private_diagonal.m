% Plot a lower bound of P^1(N_{s,t})
function res = private_diagonal(s,t,x)
     rho_B_avg = diag([s/2, (1-s-t)/2, t/2 + 1/2]);
     rho_E_avg = diag([s/2 + x/2, (1-s-t)/2 + (1-x)/2, t/2]);
     rho_E2 = diag([x, 1-x, 0]);
     res = real(VNent(rho_B_avg)-VNent(rho_E_avg)+ VNent(rho_E2)/2); 
end

% find the s,t such that Q^{(1)}(N_{s,t}) = 0.
% Under the conjecture Q^{(1)}(N_{s,t}) is optimized on the diagonal. 
function res = ci_diagonal(s,t,x)
  if (s <= 1 - s -t)
     rho_B = diag([s*x, (1-s-t)*x, t*x + 1 - x]);
     rho_E = diag([s*x, (1-s-t)*x + 1 - x, t*x]);
     res = real(VNent(rho_B)-VNent(rho_E));
  else
     rho_B = diag([s*x, (1-s-t)*x, t*x + 1 - x]);
     rho_E = diag([s*x + 1-x, (1-s-t)*x, t*x]);
     res = real(VNent(rho_B)-VNent(rho_E));
  end
end

opt = optimoptions('fminunc','disp','none');

%% 


% Define the grid for s and t
s_values = linspace(0, 1, 1000);
t_values = linspace(0, 0.5, 1000);
[S, T] = meshgrid(s_values, t_values);

% Preallocate the matrix to store private and coherent information
max_private = zeros(size(S));
max_coherent = zeros(size(S));
% Loop through the grid and find the maximum for each pair (s,t)
for i = 1:length(s_values)
    for j = 1:length(t_values)
    % Apply the additional restriction s >= 0, t >= 0, s + t <= 1, and t < 0.5
        if (S(i, j) + T(i, j) <= 1) && (T(i, j) < 0.5)
        % Define the function to be maximized for fixed s and t
            func_priv = @(x) -private_diagonal(S(i, j), T(i, j), x); % Negative for maximization
            func_ci = @(x) -ci_diagonal(S(i, j), T(i, j), x);
            % Find the maximum in the range [0, 1]
            [x_priv, f_priv] = fminbnd(func_priv, 0, 1);
            [x_ci, f_ci] = fminbnd(func_ci, 0, 1);
            max_private(i, j) = -f_priv; % Store the maximum value
            max_coherent(i, j) = -f_ci; 
        else
            max_private(i, j) = NaN;
            max_coherent(i, j) = NaN; % Set to NaN if the condition is not satisfied
        end
    end
end
%% 

% Plot the result
figure;
surf(S, T, max_private);
xlabel('s');
ylabel('t');
zlabel('P^{(1)}(N_{s,t}) and Q^{(1)}(N_{s,t})');
title('Plot of private and coherent information in the region s+t \leq 1 and t < 1/2');

alpha(0.5); % Set transparency level
grid on;
shading interp; % Smooth shading

hold on; % Retain the current plot

% Plot the second surface
mesh(S, T, max_coherent, 'FaceColor', 'none', 'EdgeColor', 'k'); % Plot of Q^1

% You can also adjust properties for better visualization
shading interp; % Apply smooth shading to the second surface as well

hold off; % Release the plot to prevent further additions