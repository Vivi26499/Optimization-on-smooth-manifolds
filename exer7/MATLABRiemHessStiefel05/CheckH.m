%% Check Hessian
% This map takes as input a symmetric matrix A, a matrix X \in St(n,p), 
% a matrix U \in T_X St(n,p) and a time step t>0. It returns the error in
% Frobenius norm between the exact hessian and the finit difference
% approximation with time step t.

function err = CheckH(A,X,U,t)

I = eye(size(A));

% Computation of Riemannian gradient at X
RG = 2 * (I-X*X')*A*X;

% Computation of the projection to T_X St(n,p) of the 
% Riemannian gradient at R_X(tU)
V = Retraction(X,t*U);
RGt = 2 * (I-V*V')*A*V;
RGtP = Proj(X,RGt);

% Finite difference approximation of the Hessian
app = (RGtP - RG) / t;

% Computation of exact hessian
ERH = RiemH(A,X,U);

% Computation of the error
err = norm(ERH - app,'fro');

end