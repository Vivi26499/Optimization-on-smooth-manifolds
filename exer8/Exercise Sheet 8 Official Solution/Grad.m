%% Computation of Riemannian gradient of f

function G = Grad(A,B,Q)
G = - 2*Q*(Q'*A*Q*B-B*Q'*A*Q);
end