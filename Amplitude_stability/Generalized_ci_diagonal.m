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

% Define the grid for s and t
s_values = linspace(0, 1, 1000);
t_values = linspace(0, 0.5, 1000);
[S, T] = meshgrid(s_values, t_values);

% Preallocate the matrix to store max values
max_values = zeros(size(S));
% Loop through the grid and find the maximum for each pair (s,t)
for i = 1:length(s_values)
    for j = 1:length(t_values)
    % Apply the additional restriction s >= 0, t >= 0, s + t <= 1, and t < 0.5
        if (S(i, j) + T(i, j) <= 1) && (T(i, j) < 0.5)
        % Define the function to be maximized for fixed s and t
            func = @(x) -ci_diagonal(S(i, j), T(i, j), x); % Negative for maximization
            % Find the maximum in the range [0, 1]
            [x_opt, fval] = fminbnd(func, 0, 1);
            max_values(i, j) = -fval; % Store the maximum value
        else
            max_values(i, j) = NaN; % Set to NaN if the condition is not satisfied
        end
    end
end
%% 

% Plot the result
figure;
surf(S, T, max_values);
xlabel('s');
ylabel('t');
zlabel('Q^{(1)}(N_{s,t})');
title('Plot of Q^{(1)}(N_{s,t}) in the region s+t \leq 1 and t < 1/2');
grid on;
shading interp; % Smooth shading
%% 

% Find the region where f(s,t) is equal to zero (within a small tolerance)
tolerance = 1e-10;
zero_region = max_values < tolerance;


%% 

% Plot the region where Q^1(N_s,t) = 0
figure;
contourf(T, S, zero_region, [1 1], 'LineColor', 'none');
xlabel('t');
ylabel('s');
title('Additivity and non-additivity properties in the region where $\mathcal{N}_{s,t}$ is neither degradable nor anti-degradable','Interpreter', 'latex', 'FontSize', 30, 'FontWeight', 'bold');
grid on;
colormap([0 0.5 0.5]); % Custom colormap for better visualization

% Add the line s + t = 1, t < 0.5 in red
hold on;
t_line = linspace(0, 0.5, 1000);
s_line = 1 - t_line;
plot(t_line, s_line, 'r', 'LineWidth', 2);
plot(t_line, zeros(size(t_line)), 'r', 'LineWidth', 2);

% Add text labels using LaTeX notation
text(0.05, 0.8, 'Additivity holds on the red boundary', 'Interpreter', 'latex', 'FontSize', 25, 'Color', 'black');
text(0.22, 0.9, '$Q^{(1)}(\mathcal{N}_{s,t} \otimes \mathcal{M}) = Q^{(1)}(\mathcal{N}_{s,t}) + Q^{(1)}(\mathcal{M})$, $\mathcal{M}$ degradable ', 'Interpreter', 'latex', 'FontSize', 30, 'Color', 'black');
text(0.22, 0.85, '$Q^{(1)}(\mathcal{N}_{s,t}) = Q(\mathcal{N}_{s,t}) $', 'Interpreter', 'latex', 'FontSize', 30, 'Color', 'black');
text(0.15, 0.5, 'Strong additivity with degradable channels fails', 'Interpreter', 'latex', 'FontSize', 30, 'Color', 'black');
text(0.15, 0.4, '$Q^{(1)}(\mathcal{N}_{s,t} \otimes \mathcal{M})>Q^{(1)}(\mathcal{N}_{s,t}) + Q^{(1)}(\mathcal{M})$', 'Interpreter', 'latex', 'FontSize', 30, 'Color', 'black');
text(0.15, 0.3, 'Log-singularity argument', 'Interpreter', 'latex', 'FontSize', 30, 'Color', 'black');
text(0.15, 0.2, 'Conjecture: Weak additivity holds', 'Interpreter', 'latex', 'FontSize', 30, 'Color', 'black');
text(0.37, 0.33, 'The region $Q^{(1)}(\mathcal{N}_{s,t}) = 0,\ P^{(1)}(\mathcal{N}_{s,t}) > 0$ ', 'Interpreter', 'latex', 'FontSize', 20, 'Color', 'black');
text(0.37, 0.3, 'Strong additivity fails', 'Interpreter', 'latex', 'FontSize', 20, 'Color', 'black');
text(0.37, 0.275, 'with higher dim erasure', 'Interpreter', 'latex', 'FontSize', 20, 'Color', 'black');
text(0.43, 0.2, 'Smith-Yard argument', 'Interpreter', 'latex', 'FontSize', 20, 'Color', 'black');
hold off;