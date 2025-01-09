function U = GradientDescent(Q, A, B, tol, max_iter)
    % GradientDescent solves the optimization problem using the Gradient Descent method
    % Inputs:
    %   Q - current point on SO(d)
    %   A, B - symmetric matrices of size d x d
    %   tol - tolerance for the stopping criterion
    %   max_iter - maximum number of iterations
    % Outputs:
    %   U - solution to the optimization problem

    b = -RiemannianGradient(Q, A, B);
    U = zeros(size(Q));
    V = b;
    for i = 1:max_iter
        H = RiemannianHessian(Q, A, B, V);
        alpha = norm(V, 'fro')^2 / trace(V' * H);
        U = U + alpha * V;
        V = V - alpha * H;
        if norm(V, 'fro') < tol * norm(b, 'fro')
            break;
        end
    end
end
