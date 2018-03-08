function T = toeplitz_adjoint(A) %#codegen
% TOEPLITZ_APPROX(A)
% For a square matrix A, find the Toeplitz approximation
% by averaging along the diagonals.
N = min(size(A));
T = complex(zeros(N,1));
T(1) = trace(A);
for n = 1:(N-1)
	tmp = diag(A,n);
    T(n+1) = 2*sum(tmp(:));
end

