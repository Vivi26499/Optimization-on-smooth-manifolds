function U = ConjugateGradient(Q, A, B, tol, max_iter)
    % ConjugateGradient solves the optimization problem using the Conjugate Gradient method
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
    W = V;
    for i = 1:max_iter
        H = RiemannianHessian(Q, A, B, W);
        alpha = norm(V, 'fro')^2 / trace(W' * H);
        disp(['Iteration ', num2str(i), ': trace(W'' * H) = ', num2str(alpha)]);
        U = U + alpha * W;
        V_old = V;
        V = V - alpha * H;
        if norm(V, 'fro') < tol * norm(b, 'fro')
            break;
        end
        beta = norm(V, 'fro')^2 / norm(V_old, 'fro')^2;
        W = V + beta * W;
    end
end