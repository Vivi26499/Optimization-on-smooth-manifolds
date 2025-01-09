%% Implementation Riemannian gradient of cost function the 
%% exercise Fréchet mean on the sphere
% Input point x on the d-1 sphere and a struct f defining f.

function gradf = gradcost(x,f)

gradf = 0;
n = f.n;
X = f.X;

for i=1:n
    d = abs(acos(x'*X(:,i)));
    gradf = gradf + (d / sin(d)) * (cos(d)*x - X(:,i));
end

gradf = gradf / n;

end