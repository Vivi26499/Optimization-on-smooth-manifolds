close all; clear; clc;
A = randn(3);
A = 0.5 * (A + A');
x0 = randn(3, 1);
x0 = x0 / norm(x0);
alpha = 1 / (2 * norm(A, "fro"));

[x, iterates] = RGDsphere(A, x0, alpha, 1000);

f = @(x) x'*A*x;
resolution = 100;
C = zeros(resolution);
[X,Y,Z] = sphere(resolution);

v = zeros(3, 1);
for i = 1:resolution
    for j = 1:resolution
        v(1) = X(i, j); 
        v(2) = Y(i, j);
        v(3) = Z(i, j); 
        C(i,j) = f(v); 
    end 
end

surf(X, Y, Z, C, "FaceAlpha", 0.9, "EdgeColor", "none");
colorbar('eastoutside');
axis equal;
hold on;
% Adding the iterates of the algorithm as red points on the sphere
plot3(iterates(1, :), iterates(2, :), iterates(3, :), 'MarkerSize', 20, 'Color', 'r');