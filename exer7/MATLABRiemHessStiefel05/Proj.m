%% Projection for Stiefel

function P = Proj(X,U)
P = U - X * (X'*U + U'*X) / 2;
end