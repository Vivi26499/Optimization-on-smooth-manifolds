%% Projection map for SO(d)

function P = Proj(Q,U)
P = Q*(Q'*U - U'*Q)/2;
end