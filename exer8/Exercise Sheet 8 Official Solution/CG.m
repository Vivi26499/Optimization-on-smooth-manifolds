%% Conjugat Gradient algorithm
% Input : g, which containes the map H and the vector b from question 7, 
% a tolerence for the convergence of CG and a maximum number of iterations.
% Outup : Approximation of minimizer of the map g from question 7.
% ATTENTION : H(U) and b have to be matrices or vectors in R^n
% H has to be s.p.d.

function V = CG(g,tol,iter)

% Initialization
H = @(U) g.H(U);
b = g.b;
V = zeros(size(b));
Rl = b;
Rn = b;
P = Rn;

for i=1:iter
    Rl = Rn;
    Hint = H(P);
    alpha = (Rl(:)'*Rl(:)) / (P(:)'*Hint(:));
    V = V + alpha * P;
    Rn = Rl - alpha * Hint;

    if norm(Rn) < tol
        fprintf('Algorithm did converge\n')
        break;
    end

    beta = (Rn(:)'*Rn(:)) / (Rl(:)'*Rl(:));
    P = Rn + beta * P;

end

fprintf('Norm of the residual\n')
norm(Rn)

end
