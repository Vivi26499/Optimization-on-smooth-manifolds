function diff = checkH(A, X, U, t, retraction_type)
    % Define sym and Proj functions
    sym = @(A) 0.5 * (A + A');
    Proj = @(X, U) U - X * sym(X' * U);

    % Calculate the Riemannian Hessian of f
    function H = riemannian_hessian(A, X, U)
        H = 2 * Proj(X, A * U) - 2 * Proj(X, U * (X' * A * X));
    end

    % Calculate the Riemannian gradient of f
    function g = riemannian_gradient(A, X)
        g = 2 * (eye(size(X, 1)) - X * X') * A * X;
    end

    % Define the retraction functions
    function RX = qr_retraction(X, V)
        [Q, R] = qr(X + V, 0);
        RX = Q * diag(sign(diag(R)));
    end

    function RX = metric_projection(X, V)
        [U, ~, V] = svd(X + V, 0);
        RX = U * V';
    end

    function RX = polar_retraction(X, V)
        RX = (X + V) * (eye(size(V, 2)) + V' * V)^(-1/2);
    end

    % Select the retraction method
    switch retraction_type
        case 'qr'
            retraction = @qr_retraction;
        case 'metric'
            retraction = @metric_projection;
        case 'polar'
            retraction = @polar_retraction;
        otherwise
            error('Unknown retraction type');
    end

    % Calculate the finite difference approximation of the Riemannian Hessian
    function H_approx = finite_difference_hessian(A, X, U, t)
        RX = retraction(X, t * U);
        grad_f_RX = riemannian_gradient(A, RX);
        grad_f_X = riemannian_gradient(A, X);
        H_approx = (Proj(X, grad_f_RX) - grad_f_X) / t;
    end

    % Compute the Riemannian Hessian and the finite difference approximation
    H = riemannian_hessian(A, X, U);
    H_approx = finite_difference_hessian(A, X, U, t);

    % Calculate the difference
    diff = norm(H - H_approx, 'fro');
end