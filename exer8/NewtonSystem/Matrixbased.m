function U = Matrixbased(Q, A, B)
    % Matrixbased solves the Newton system for optimization on SO(d)
    % Inputs:
    %   Q - current point on SO(d)
    %   A, B - symmetric matrices of size d x d
    % Outputs:
    %   U - solution to the Newton system

    d = size(Q, 1);
    n = d * (d - 1) / 2;  % dimension of SO(d)
    
    % Generate an orthonormal basis for T_Q SO(d)
    Z = OrthonormalBasis(d, Q);

    % The Newton system
    H = zeros(n, n);
    b = zeros(n, 1);
    for i = 1:n
        for j = 1:n
            H(i, j) = trace(Z(:, :, i)' * RiemannianHessian(Q, A, B, Z(:, :, j)));
        end
        b(i) = -trace(Z(:, :, i)' * RiemannianGradient(Q, A, B));
    end
    
    % Solve the Newton system
    alpha = H \ b;
    U = sum(alpha' .* reshape(Z, d*d, n), 2);
    U = reshape(U, d, d);
end

function Z = OrthonormalBasis(d, Q)
    % OrthonormalBasis generates an orthonormal basis for T_Q SO(d)
    % Inputs:
    %   d - dimension of the space
    %   Q - current point on SO(d)
    % Output:
    %   Z - 3D array containing the orthonormal basis matrices

    Z = zeros(d, d, d * (d - 1) / 2);
    idx = 1;
    for i = 1:d-1
        for j = i+1:d
            E = zeros(d, d);
            E(i, j) = 1;
            E(j, i) = -1;
            Z(:, :, idx) = Q * E / sqrt(2);
            idx = idx + 1;
        end
    end
end