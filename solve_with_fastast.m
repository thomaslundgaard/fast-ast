function [x, u, v, info] = solve_with_fastast(y, tau, varargin)
%SOLVE_WITH_FASTAST  Fast Atomic Norm Soft Thresholding.
%   output = solve_with_fastast(y, tau, ...) solves the atomic norm soft
%   thresholding problem based on the method presented in [1].
%
%	[1] T. L. Hansen and T. L. Jensen - "A Fast Interior Point Method for
%	Atomic Norm Soft Thresholding", submitted to IEEE Transactions on Signal
%	Processing, 2018
%   
%   y is vector with time domain measurements.
%
%	tau is the tradeoff parameter.
%
%   Further arguments can be specified as name,value pairs, e.g. as:
%     solve_with_fastast(y,index,N, 'verbose', true, ...)
%   Each entry in the 'as' structure specified in the source code can be given this
%   way. If nothing is specified, the default value specified in the source
%   code is used.
%
%	Copyright (c) 2018
%	Contributors:
%		Thomas Lundgaard Hansen
%		Tobias Lindstrøm Jensen
%
%	This software is being released under the MIT license (see LICENSE file).

	N = length(y);

	% Add path if needed
	if exist('fastast_worker.m')~=2
		mfilepath=fileparts(which(mfilename));
		addpath(fullfile(mfilepath, 'util'));
	end

	% First define everyting that can be set and parse input to get method
	as = struct();
	as.verbose = false;
	as.test = false;
	as.method = 1;			% 1 => L-BFGS
							% 2 => Newton's method (non-optimized O(N^4) implementation)
	as.tol_abs = 1e-4;		% Stop when duality gap is below this value
	as.tol_rel = 1e-4;		% or when relative duality gap is below this value
	as.max_iters = 2000;
	as.mu = 2;				% How much to update t
	as.stepsize_reduction = 1/4;
	as.num_saved_lbfgs = 2*N-1;
	as.w = [2; zeros(2*N-2,1)];
	as.u0 = [10*var(y,1); zeros(2*N-2, 1)];
	as = parse_varargin(as, varargin);

	% If method is 2, use another set of defaults
	if as.method == 2
		as.tol_abs = 1e-7;	% Stop when duality gap is below this value
		as.tol_rel = 1e-7;	% Or when relative duality gap is below this value
		as.max_iters = 200;
		as.mu = 10;
		as.stepsize_reduction = 7/10;

		% Parse again to allow these defaults to be overwritten
		as = parse_varargin(as, varargin);
	end

	% Call fastast_solver.
	% We typecast some variables to adhere to the prototype used in the mex version.
	[x, u, v, info] = fastast_worker(y, tau, as.w, as.u0, logical(as.verbose), int32(as.method), ...
							as.tol_abs, as.tol_rel, int32(as.max_iters), as.mu, ...
							as.stepsize_reduction, int32(as.num_saved_lbfgs));
end

