function [x, iterates] = RGDsphere(A, x0, alpha, maxiter)
conv = false;
x = x0;
iterates = zeros(length(x), maxiter);
iterates(:, 1) = x;
epsilon = 1e-15;

for k = 1:(maxiter - 1)
    grad = 2 * (eye(length(x)) - x * x') * A * x;
    v = -alpha * grad;
    x = (x + v) / norm(x + v);
    iterates(:, k+1) = x;
    
    if norm(grad) <= epsilon
        conv = true;
        iterates = iterates(:, 1:k+1);
        disp('RGD converges');
        disp('Final position: ')
        disp(x);
        disp(['Final value: ', num2str(x' * A * x)]);
        disp(['Number of iterations: ', num2str(k)]);
        break;
    end
end

if conv == false
    disp('RGD did not converge');
    disp('Final position: ')
    disp(x);
    disp(['Final value: ', num2str(x' * A * x)]);
    disp(['Number of iterations: ', num2str(maxiter)]);
end

end