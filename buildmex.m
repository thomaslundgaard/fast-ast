function buildmex(force_no_recursion)
	% Set largest problem size (i.e. length(y)) which will be supported by the generated mex.
	N = 2^12;
	cd util

	if nargin<1 || ~islogical(force_no_recursion)
		 force_no_recursion = false;
	end

	if ~verLessThan('matlab', '9.0') && ~force_no_recursion
		% This builds the complete algorithms as mex for the fastest
		% performance.
		fprintf('MATLAB version >= 9.0, so building the fastest version.\n');
		fprintf('(Based on recursion in mex version of genschur).\n');

		% Build FastAST
		fprintf('Building fastast_worker...\n');
		args = {};
		args{end+1} = coder.typeof(complex(1,1), [N 1], [true false]);	% y
		args{end+1} = coder.typeof(1);									% tau
		args{end+1} = coder.typeof(1, [2*N-1 1], [true false]);			% w
		args{end+1} = coder.typeof(1, [2*N-1 1], [true false]);			% u0
		args{end+1} = coder.typeof(true);								% verbose
		args{end+1} = coder.typeof(int32(1));							% method
		args{end+1} = coder.typeof(1);									% tol_abs
		args{end+1} = coder.typeof(1);									% tol_rel
		args{end+1} = coder.typeof(int32(1));							% max_iters
		args{end+1} = coder.typeof(1);									% mu
		args{end+1} = coder.typeof(1);									% stepsize_redution
		args{end+1} = coder.typeof(int32(1));							% num_saved_lbfgs
		codegen -o fastast_worker fastast_worker.m -args args

	else
		% This builds only the computation-heavy subfunctions as mex. The main
		% algorithms are still run in matlab. This may in fact be faster than
		% building everything into mex if there is multiple threads available
		% for matlab. It is also the only thing which works on older vesion of
		% matlab.
		fprintf('MATLAB version < 9.0, so building the slower version.\n');
		fprintf('(No recursion in mex version of genschur).\n');

		% Uses recursion, so only works in never versions of matlab
		%fprintf('Building genschur...\n');
		%args = {};
		%args{end+1} = coder.typeof(complex(1,1), [1 N], [false true]);
		%args{end+1} = coder.typeof(complex(1,1), [1 N], [false true]);
		%codegen -o genschur genschur.m -args args

		% As a fallback we can use mex only for the innermost iteration of
		% the generalized Schur algorithm. This will be slower than the
		% above.
		Nrecur = 32;
		fprintf('Building genschur_recurrence...\n');
		args = {};
		args{1} = coder.typeof(complex(1,1), [1 Nrecur], [false true]);
		args{2} = coder.typeof(complex(1,1), [1 Nrecur], [false true]);
		codegen -o genschur_recurrence genschur_recurrence.m -args args

		% Build levin
		fprintf('Building levin...\n');
		args = {};
		args{end+1} = coder.typeof(complex(1,1), [N 1], [true false]);
		codegen -o levin levin.m -args args

		% Build trTinvTinv
		fprintf('Building trTinvTinv...\n');
		args = {};
		args{end+1} = coder.typeof(complex(1,1), [N 1], [true false]);
		args{end+1} = coder.typeof(1);
		codegen -o trTinvTinv trTinvTinv.m -args args

		% The un-mexed version of these is faster (on my MacBook at least)
		%% Build form_barrier_hessian
		%fprintf('Building form_barrier_hessian...\n');
		%args = {};
		%args{end+1} = coder.typeof(complex(1,1), [N N], [true true]);
		%args{end+1} = coder.typeof(complex(1,1), [N 1], [true false]);
		%codegen -o form_barrier_hessian form_barrier_hessian.m -args args

		%% Build lbfgs_get_direction
		%fprintf('Building lbfgs_get_direction...\n');
		%% Gives segfault on old version of matlab for some reason. On old
		%% version we just use the raw matlab verion.
		%args = {};
		%args{end+1} = struct();
		%args{end}.max_saved = coder.typeof(int32(1));
		%args{end}.num_saved = coder.typeof(int32(1));
		%args{end}.newest = coder.typeof(int32(1));
		%args{end}.ss = coder.typeof(1, [2*N-1 2*N-1], [true true]);
		%args{end}.yy1 = coder.typeof(1, [2*N-1 2*N-1], [true true]);
		%args{end}.yy2 = coder.typeof(1, [2*N-1 2*N-1], [true true]);
		%args{end}.sy1 = coder.typeof(1, [2*N-1 1], [true false]);
		%args{end}.sy2 = coder.typeof(1, [2*N-1 1], [true false]);
		%args{end+1} = coder.typeof(1, [2*N-1 1], [true false]);
		%args{end+1} = coder.typeof(1, [2*N-1 1], [true false]);
		%args{end+1} = coder.typeof(1);
		%codegen -o lbfgs_get_direction lbfgs_get_direction.m -args args
	end

	% Build admm_worker
	% The un-mexed version is faster (on my MacBook at least)
	%fprintf('Building admm_worker...\n');
	%args = {};
	%args{end+1} = coder.typeof(complex(1,1), [N 1], [true false]);	% y
	%args{end+1} = coder.typeof(1);									% tau
	%args{end+1} = coder.typeof(true);								% verbose
	%args{end+1} = coder.typeof(int32(1));							% maxIter
	%args{end+1} = coder.typeof(1);									% tol_abs
	%args{end+1} = coder.typeof(1);									% tol_rel
	%codegen -o admm_worker admm_worker.m -args args

	% Build toeplitz_adjoint (only used by admm_worker)
	fprintf('Building toeplitz_adjoint...\n');
	args = {};
	args{end+1} = coder.typeof(complex(1,1), [N N], [true true]);
	codegen -o toeplitz_adjoint toeplitz_adjoint.m -args args

	cd ..
end

