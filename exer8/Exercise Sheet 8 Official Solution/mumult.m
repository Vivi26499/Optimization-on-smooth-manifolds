%% Implementation of the mu-mode product
% Should just be used for small tesnors !
% Input : tensor A, matrix B and mode mu
% Output : tensor C

function C = mumult(A,B,mu)

sA = size(A);
sB = size(B);
s_mu = sA(1,mu);

% check sizes
if s_mu ~= sB(1,2)

    fprintf('there is a dimensional problem \n')
    return

else


A_mu = ten2mum(A,mu); % mu-mode matricization of the tensor
C_mu = B * A_mu; % matrix product with mu-mode matricization
C = mum2ten(C_mu,sA,mu); % mu-mode tensorization of C_mu

end

end