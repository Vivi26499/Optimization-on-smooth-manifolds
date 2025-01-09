%% Gradient Descent algorithm
% Input : g, which containes the map H and the vector b from question 7, 
% a tolerence for the convergence of GD and a maximum number of iterations.
% Outup : Approximation of minimizer of the map g from question 7.
% ATTENTION : H(U) and b have to be matrices or vectors in R^n
% H has to be s.p.d.

function V = GD(g,tol,iter)

% Initialization
H = @(U) g.H(U);
b = g.b;
V = zeros(size(b));
R = b;

for i=1:iter
    Hint = H(R);
    alpha = R(:)'*R(:) / (R(:)'*Hint(:));
    V = V + alpha * R;
    R = R - alpha * Hint;

    if norm(R) < tol
        fprintf('Algorithm did converge\n')
        break;
    end

end

fprintf('Norm of the residual\n')
norm(R)

end

