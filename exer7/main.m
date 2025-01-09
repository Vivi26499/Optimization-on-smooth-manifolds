% Define the matrix A and the point X on the Stiefel manifold
n = 5; % Example dimension
p = 3; % Example dimension
A = randn(n, n);
A = 0.5 * (A + A'); % Ensure A is symmetric
X = orth(randn(n, p)); % Ensure X is on the Stiefel manifold

% Define the tangent vector U
U = randn(n, p);

% Define the step sizes and retraction types
t_values = [10^-1, 10^-2, 10^-4, 10^-8];
retraction_types = {'qr', 'metric', 'polar'};

% Initialize the results table
results = zeros(length(t_values), length(retraction_types));

% Loop over the step sizes and retraction types and compute the difference
for i = 1:length(t_values)
    t = t_values(i);
    for j = 1:length(retraction_types)
        retraction_type = retraction_types{j};
        diff = checkH(A, X, U, t, retraction_type);
        results(i, j) = diff;
        fprintf('Difference for t = %e and retraction = %s: %f\n', t, retraction_type, diff);
    end
end

% Display the results table
disp('Results Table:');
disp(array2table(results, 'VariableNames', retraction_types, 'RowNames', arrayfun(@(t) sprintf('t=%e', t), t_values, 'UniformOutput', false)));