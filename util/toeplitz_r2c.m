function out = toeplitz_r2c(u)
	N = (length(u)+1) / 2;
	u(1) = 2*u(1);
	out = toeplitz( u(1:N) + [0;1j*u(N+1:end)] );
end

