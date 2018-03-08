% Copyright (c) 2018: Thomas Lundgaard Hansen, Tobias Lindstrøm Jensen
% This software is being released under the MIT license (see LICENSE file). 
function [x, u, v, info] = solve_with_cvx(y, tau, varargin)
	N = length(y);

	% Add path if needed
	if exist('fastanm_worker.m')~=2
		mfilepath=fileparts(which(mfilename));
		addpath(fullfile(mfilepath, 'util'));
	end

	% First define everyting that can be set and parse input to get method
	as = struct();
	as.verbose = false;
	as.precision = 'default';
	as.solver = 'SeDuMi';
	as.w = [2; zeros(2*N-2,1)];
	as = parse_varargin(as, varargin);

	cvx_begin
		cvx_solver(as.solver);
		if as.verbose
			cvx_quiet false
		else
			cvx_quiet true
		end
		cvx_precision(as.precision);

		variable T(N, N) hermitian toeplitz
		variable v(1)
		variable x(N) complex
		variable e(N) complex
		dual variable z

		minimize pow_pos(norm(e), 2) ...
			+ tau*v ...
			+ tau*T(1,1)*as.w(1)/2 ...
			+ tau*real(T(1,2:end))*as.w(2:N) ...
			+ tau*imag(T(1,2:end))*as.w(N+1:end);
		subject to
			[T x; x' v] == hermitian_semidefinite(N+1)
			z : e == x - y
	cvx_end

	info = cvx_status;
	u = T(1,:).';
	u = [real(u(1))/2; real(u(2:end)); imag(u(2:end))];
end

