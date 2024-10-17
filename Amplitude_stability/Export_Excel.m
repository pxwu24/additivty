

% Assume x, y, z are your high-dimensional vectors
% Combine the vectors into a matrix where each vector is a column
% Ensure the vectors are column vectors
x = gamma(:);  % Converts x to a column vector
y = eta(:);  % Converts y to a column vector
z = ratio_array(:);  % Converts z to a column vector


data = [x,y,z];

% Write the matrix to an Excel file
writematrix(data, 'output.xlsx');

% Optionally, specify the sheet and range
writematrix(data, 'output.xlsx', 'Sheet', 'Sheet1', 'Range', 'A1');