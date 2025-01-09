%% Implementation of the mu-mode matricization A to A_\mu
% Should just be used for small tesnors !
% Input : tensor A and mode mu
% Output : mu-mode of A

function A_mu = ten2mum(A,mu)

% Size of the tensor
s = size(A);
d = max(size(s));

% Control if the entry is acceptable
    if mu > d

        fprintf('there is a dimensional problem \n')
        return

    end

% Control of trivial case d = 1
if d == 1

    A_mu = A;

else
    
    n_mu = s(1,mu); % size of the mu-mode
    r_mu = prod(s) / n_mu; % size of the columns of A_mu
    p_mu = 1:d;
    p_mu(1,1) = mu;
    p_mu(1,mu) = 1; % permuted indices

    A = permute(A,p_mu); % permuted tensor
    A_mu = reshape(A,n_mu,r_mu); % 1-mode matricization

end

end