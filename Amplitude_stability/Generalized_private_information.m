dim = [a1,b1,e1,3];
% number of partitions
resolution = 100;
% Create arrays to store valid (s, t) pairs and corresponding Q^(1)(N_{s,t})
s_values = [];
t_values = [];
p_values = [];
temp_private = 0;
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
            temp_private = optimize_private_information(V,dim,3);
            disp(['The parameters of generalized Platypus channel are: ', 's= ', num2str(s),', t= ', num2str(t)])
            disp(['private information of generalized Platypus channel = ', num2str(temp_private)])
            p_values = [p_values; temp_private]; % Store the function value P^(1)(N_{s,t})
        end
    end
end

% Create a scatter plot of the points in the region
scatter3(s_values, t_values, p_values, 36, p_values, 'filled');

% Label axes
xlabel('s');
ylabel('t');
zlabel('P^(1)(N_{s,t})');
title('Plot of P^(1)(N_{s,t}) in the region s+t \leq 1 and t < 1/2');

% Add a color bar to indicate the function value
colorbar;

% Set view for better visualization
view(3);


%% % Ensure the vectors are column vectors
x = s_values(:);  % Converts x to a column vector
y = t_values(:);  % Converts y to a column vector
z = p_values(:);  % Converts z to a column vector


data = [x,y,z];

% Write the matrix to an Excel file
writematrix(data, 'private_info.xlsx');

%% 
