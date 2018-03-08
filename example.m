clear all
addpath util


if 0
	% To replay a run, first save it with
	% save('saved.mat', 'xs','y','N','M','o','c','sigma','w');
	% Then we can load it here
	load saved.mat
else
	%rng(2);

	% Parameters
	snr = 20;
	N = 2^6;
	M = max(1,round(N/10));

	% Generate signal
	n = (0:N-1)';
	%o = [0.1 0.11 0.52];
	o = sort(rand(M,1));
	c = randn(length(o), 2) * [1;1j];
	c = c./abs(c);

	xs = zeros(N, 1);
	for m = 1:length(o)
		xs = xs + c(m)*exp(1j*2*pi*o(m)*n);
	end

	noiseVar = mean(abs(xs).^2) / 10^(snr/10);
	w = sqrt(noiseVar/2) * randn(N, 2) * [1;1j];
	y = xs + w;
end

% Prepare stuff
tau = sqrt(noiseVar)*(1+1/log(N))*sqrt(N*log(N) + N*log(4*pi*log(N)));
f = @(v,x,u) norm(x-y)^2 + tau*v + tau * 2 * real(u(1));

% Solve
tic
[xx, uu, vv, info] = solve_with_fastast(y, tau, 'verbose',1, 'method',1);
toc

%tic
%[xxx, uuu, vvv, info] = solve_with_admm(y, tau);
%toc

tic
[x, u, v, info] = solve_with_cvx(y, tau);
toc

% Recover frequencies
%omega = recover_via_esprit(toeplitz_r2c(u), length(o));
omega = recover_via_dual(x, y, tau);
%omega2 = recover_via_esprit(toeplitz_r2c(uu), length(o));
omega2 = recover_via_dual(xx, y, tau);

% Print results
fprintf('FastAST\tCVX\tDifference\n');
if length(omega2) == length(omega)
	omega_ = [omega2, omega, abs(omega2-omega)]
else
	omega
	omega2
end
obj_ = [f(vv,xx,uu), f(v,x,u), f(vv,xx,uu) - f(v,x,u)]
max_u_uu = max(abs(u-uu))

