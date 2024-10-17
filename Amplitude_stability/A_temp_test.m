% Define the channel parameters
V = rand(8, 8);  % Random example matrix of appropriate size
dim = [2, 2, 2, 2];  % Example dimensions for subsystems
k = 3;  % Number of Kraus operators

% Run the function with the default optimization settings
[pi, p, R] = optimize_private_information(V, dim, k);
