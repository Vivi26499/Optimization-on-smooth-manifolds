%% Computation of Riemannian Hessian of f

function H = Hess(A,B,Q,U)
Hint = 2*(Q*B*Q'*A*U+Q*B*U'*A*Q+U*B*Q'*A*Q - A*U*B);
H = Proj(Q,Hint);
end