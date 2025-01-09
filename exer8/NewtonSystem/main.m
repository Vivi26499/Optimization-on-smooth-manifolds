% main.m
% Compare the performance of ConjugateGradient, GradientDescent, and Matrixbased methods
% for different dimensions of Q

% Parameters
tol = 1e-6;
max_iter = 100;

% Dimensions to test
dimensions = [3, 10];

% Initialize results
results = struct('dimension', [], 'ConjugateGradient', [], 'GradientDescent', [], 'Matrixbased', [], 'CG_Error', [], 'GD_Error', []);

for d = dimensions
    % Generate random symmetric matrices A and B
    A = randn(d);
    A = (A + A') / 2;
    B = randn(d);
    B = (B + B') / 2;
    
    % Generate a random orthogonal matrix Q and ensure it is a rotation matrix
    [Q, ~] = qr(randn(d));
    
    % Measure time for ConjugateGradient
    tic;
    U_CG = ConjugateGradient(Q, A, B, tol, max_iter);
    time_CG = toc;
    
    % Measure time for GradientDescent
    tic;
    U_GD = GradientDescent(Q, A, B, tol, max_iter);
    time_GD = toc;
    
    % Measure time for Matrixbased
    tic;
    U_MB = Matrixbased(Q, A, B);
    time_MB = toc;
    
    % Calculate errors
    error_CG = norm(U_CG - U_MB, 'fro');
    error_GD = norm(U_GD - U_MB, 'fro');
    
    % Store results
    results(end+1) = struct('dimension', d, 'ConjugateGradient', time_CG, 'GradientDescent', time_GD, 'Matrixbased', time_MB, 'CG_Error', error_CG, 'GD_Error', error_GD);
end

% Display time results
disp('Time Results:');
disp('Dimension | ConjugateGradient | GradientDescent | Matrixbased');
for i = 1:length(results)
    fprintf('%9d | %17.6f | %15.6f | %11.6f\n', results(i).dimension, results(i).ConjugateGradient, results(i).GradientDescent, results(i).Matrixbased);
end

% Display error results
disp('Error Results:');
disp('Dimension | CG_Error | GD_Error');
for i = 1:length(results)
    fprintf('%9d | %8.6e | %8.6e\n', results(i).dimension, results(i).CG_Error, results(i).GD_Error);
end