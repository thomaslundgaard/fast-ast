% Copyright (c) 2018: Thomas Lundgaard Hansen.
% This software is being released under the MIT license (see LICENSE file). 
function H = form_barrier_hessian(R, delta)
	%#codegen
	N = length(delta);
	Nfft = 2^ceil(log2(2*N - 1));
	assert(N<=4096);	% allow static memory allocation
	assert(Nfft<=8192);	% allow static memory allocation
	R = bsxfun(@times, R, 1./sqrt(delta(:)'));
	S = fft(R, Nfft);
	M = complex(zeros(Nfft,Nfft));

	% Carry out sums
	%for i = 1:N
		%M = M + S(:,i) * S(:,i)';
	%end
	M = S*S';

	% Form whats multiplied by W matrices
	M = M .* M.';

	% Form conj(W)*M
	M = fft(M);
	M = M(1:N,:)/Nfft;

	% Form A
	A = ifft(M')';
	A = 2*A(:,1:N);

	% Form B
	B = ifft(M.').';
	B = 2*B(:,1:N);

	% Form Hessian
	H = [real(A+B), ...
			real(-1j*A(:,2:end) + 1j*B(:,2:end)); ...
			real(-1j*A(2:end,:) - 1j*B(2:end,:)), ...
			real(-A(2:end,2:end) + B(2:end,2:end))];
end

