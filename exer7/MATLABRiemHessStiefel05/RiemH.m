%% Computation of Riemannian Hessian for the map 
%% X \in St(n,p) \mapsto Tr(X'AX) in direction U

function RH = RiemH(A,X,U)
I = eye(size(X*X'));
RH = 2*A*U - 2*X*Sym(X'*A*U) - (I - X*X' / 2) * (U*X'*A*X + U*X'*A*X) + X*X'*A*X*U'*X;
end