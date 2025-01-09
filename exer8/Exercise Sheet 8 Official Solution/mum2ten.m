%% Implementation of the mu-mode matricization A_\mu to A
% Should just be used for small tesnors !
% Input : mu-mode  of a tensor A, dimension of tensor (index) and mode mu
% Output : A

function A = mum2ten(A_mu,index,mu)

% Size of the tensor
s = size(A_mu);
d = max(index);

% Control if the entry is acceptable
if (mu > d) || (prod(index) / index(1,mu) ~= s(1,2))

   fprintf('there is a dimensional problem \n')
   return

end

% Control of trivial case d = 1
if d == 1

    A = A_mu;

else

    n_mu = index(1,mu);
    p_mu = 1:d;
    p_mu(1,1) = mu;
    p_mu(1,mu) = 1; % permuted indices
    p_index = index;
    p_index(1,mu) = index(1,1);
    p_index(1,1) = n_mu; % permuted index set

    A = reshape(A_mu,p_index); % permuted tensor
    A = permute(A,p_mu); % non-permuted tensor

end

end