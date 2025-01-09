function riem_grad_f = RiemannianGradient(Q, A, B)
    % RiemannianGradient computes the Riemannian gradient of f at Q
    % Inputs:
    %   Q - current point on SO(d)
    %   A, B - symmetric matrices of size d x d
    % Output:
    %   riem_grad_f - Riemannian gradient at Q

    euc_grad_f = 2 * (Q * B * Q' * A * Q - A * Q * B);
    riem_grad_f = ProjectToTangentSpace(Q, euc_grad_f);
end