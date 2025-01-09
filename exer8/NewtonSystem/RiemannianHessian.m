function riem_hess_f_Z = RiemannianHessian(Q, A, B, Z)
    % RiemannianHessian computes the Riemannian Hessian of f at Q along Z
    % Inputs:
    %   Q - current point on SO(d)
    %   A, B - symmetric matrices of size d x d
    %   Z - direction in the tangent space T_Q SO(d)
    % Output:
    %   riem_hess_f_Z - Riemannian Hessian at Q along Z

    euc_hess_f_Z = 2 * (Z * B * Q' * A * Q + Q * B * Z' * A * Q + Q * B * Q' * A * Z - A * Z * B);
    riem_hess_f_Z = ProjectToTangentSpace(Q, euc_hess_f_Z);
end