function out = toeplitz_r2c_adjoint(X);
%
% The adjoint of a Toeplitz matrix
% The input is a complex NxN matrix.
% The output is a real 2N-1 vector
%

N = size(X, 1);

u = complex(zeros(N, 1));

for k=1:N
    u(k) = 2*sum(diag(X, k-1));
end

out = [real(u); imag(u(2:end))];

