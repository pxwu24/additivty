disp(length(s_values))
disp(length(t_values))
disp(length(f_values))

s_temp = s_values;
t_temp = t_values;
f_temp = f_values;

% Create a grid for interpolation
numPoints = 100;  % Number of points to interpolate over
xlin = linspace(min(s_temp), max(s_temp), numPoints);
ylin = linspace(min(t_temp), max(t_temp), 50);
[X, Y] = meshgrid(xlin, ylin);

% Interpolate the scattered data to the grid
Z = griddata(s_temp, t_temp, f_temp, X, Y, 'cubic');  % 'cubic' interpolation for smooth surface

% Plot the surface
figure;
surf(X, Y, Z);
xlabel('s');
ylabel('t');
zlabel('Q^{(1)}(N_{s,t})');
title('Plot of $Q^{(1)}(\mathcal{N}_{s,t})$ for $s+t \le 1$ and $t < 1/2$', 'Interpreter', 'latex');
