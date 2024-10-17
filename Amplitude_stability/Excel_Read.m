% Read the entire sheet from the Excel file
data = readmatrix('private_info.xlsx');

% Extract the first three columns as vectors
vector1 = data(:, 1);
vector2 = data(:, 2);
vector3 = data(:, 3);
%% 
value = [];
for i = 1:length(vector1)
    if (vector1(i) == 0.3) && (vector2(i) == 0.4)
        value = [value,vector3(i)];
    end
end
disp(value)

%% 

% Create a scatter plot of the points in the region
scatter3(vector1, vector2, vector3, 36, vector3, 'filled');

% Label axes
xlabel('s');
ylabel('t');
zlabel('P^{(1)}(N_{s,t})');
title('Plot of P^{(1)}(N_{s,t}) in the region s+t \leq 1 and t < 1/2');

% Add a color bar to indicate the function value
colorbar;

% Set view for better visualization
view(3);

%% % Convert x, y, z to matrices suitable for surf
[Xq, Yq] = meshgrid(unique(vector1), unique(vector2));
Zq = griddata(vector1, vector2, vector3, Xq, Yq);
% Plot the resulting meshgrid as a surface plot
surf(Xq, Yq, Zq);


xlabel('s');
ylabel('t');
zlabel('P^{(1)}(N_{s,t})');
title('Plot of P^{(1)}(N_{s,t}) in the region s+t \leq 1 and t < 1/2');

% Add a color bar to indicate the function value
colorbar;

% Set view for better visualization
view(3);