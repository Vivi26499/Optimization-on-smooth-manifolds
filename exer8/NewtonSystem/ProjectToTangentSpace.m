function P = ProjectToTangentSpace(Q, U)
    % ProjectToTangentSpace projects U onto the tangent space T_Q SO(d)
    % Inputs:
    %   Q - current point on SO(d)
    %   U - matrix to be projected
    % Output:
    %   P - projection of U onto T_Q SO(d)

    P = Q * skew(Q' * U);
end

function S = skew(M)
    % skew computes the skew-symmetric part of a matrix
    % Inputs:
    %   M - input matrix
    % Output:
    %   S - skew-symmetric part of M

    S = (M - M') / 2;
end