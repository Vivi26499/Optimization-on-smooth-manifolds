%% Function sym which gives the symetric version of a matrix

function S = Sym(A)
S = (A + A') / 2;
end