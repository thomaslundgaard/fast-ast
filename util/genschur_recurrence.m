%% Generalized Schur Algorithm
%% Uses the recurrence relations in the beginning of [1, Sec 3.1].
%% Takes (n-1)'th degree polynomials alpha_0 and beta_0, such that
%% phi_0=alpha_0/beta_0 is a Schur function.
%% Returns the n'th Schur polynomials (of degree n-1) xi_n, eta_n, obtained by
%% running n steps of the Schur algorithm on phi_0.
%% Based on [1] G. S. Ammar and W. B. Gragg, "The Generalized Schur Algorithm
%% for the Superfast Solution of Toeplitz Systems"
% Copyright (c) 2018: Thomas Lundgaard Hansen.
% This software is being released under the MIT license (see LICENSE file). 
function [xi_n, eta_n, gamma] = genschur_recurrence(alpha_0, beta_0) %#codegen
	n = length(alpha_0);
	assert(n<=32);	% allow static memory allocation

	gamma = complex(ones(n,1), zeros(n,1));
	xi_n = complex(zeros(1,n), zeros(1,n));
	eta_n = complex(zeros(1,n), zeros(1,n));
	alpha = alpha_0;
	beta = beta_0;
	for k = 1:n
		gamma(k) = alpha(1) / beta(1);
		if 1-abs(gamma(k))^2 <= 0
			% This matrix is not positive definite
			return;
		end
		alpha_new = alpha(2:n-k+1) - gamma(k) * beta(2:n-k+1);
		beta(1:n-k) = beta(1:n-k) - gamma(k)' * alpha(1:n-k);
		alpha(1:n-k) = alpha_new;
		if k == 1
			xi_n(1) = gamma(k);
			eta_n(1) = 1;
		else
			xi_new = [0, gamma(k)*conj(eta_n(k-1:-1:1))] ...
				+ [xi_n(1:k-1), 0];
			eta_n(1:k) = [0, gamma(k)*conj(xi_n(k-1:-1:1))] ...
				+ [eta_n(1:k-1), 0];
			xi_n(1:k) = xi_new;
		end
	end
end

