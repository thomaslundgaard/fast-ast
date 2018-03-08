function [x, u, v, info] = admm_solver(y, tau, varargin)
	N = length(y);

	if exist('admm_worker.m')~=2
		mfilepath=fileparts(which(mfilename));
		addpath(fullfile(mfilepath, 'util'));
	end
	as = struct();
	as.verbose = false;
	as.max_iters = 5000;	% maximum number of ADMM steps
	as.tol_abs = 1e-4;		% absolute tolerance 1e-3
	as.tol_rel = 1e-5;		% iteration level relative tolerance 1e-4
	as = parse_varargin(as, varargin);

	[x, u, v, info] = admm_worker(y, tau, ...
		logical(as.verbose), int32(as.max_iters), ...
		as.tol_abs, as.tol_rel);
end

