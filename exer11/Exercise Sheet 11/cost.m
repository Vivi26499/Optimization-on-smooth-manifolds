%% Implementation cost function for the exercise Fréchet mean on the sphere
% Input point x on the d-1 sphere and a struct f defining f.

function fval = cost(x,f)

fval = 0;
n = f.n;
X = f.X;

for i=1:n
    d2 = acos(x'*X(:,i))^2;
    fval = fval + d2;
end

fval = fval / (2*n);

end