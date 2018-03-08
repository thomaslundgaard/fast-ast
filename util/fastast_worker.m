% Copyright (c) 2018: Thomas Lundgaard Hansen, Tobias Lindstrøm Jensen.
% This software is being released under the MIT license (see LICENSE file). 
function [x, u, v, info] = fastast_worker(y, tau, w, u0, verbose, method, ...
							tol_abs, tol_rel, max_iters, mu, ...
							stepsize_reduction, num_saved_lbfgs)
	
	% Print warnings regarding convergence, numerical issues, etc?
	print_warnings = true;

	% Predefine possible output values for codegen
	N = length(y);
	info = struct();
	info.iters = -1;
	info.gap = nan(1, max_iters);
	info.relgap = nan(1, max_iters);
	info.objective = nan(1, max_iters);
	info.g = nan(1, max_iters);
	info.t = nan(1, max_iters);
	info.v = nan(1, max_iters);
	info.x = complex(nan(N, max_iters));
	info.u = nan(2*N-1, max_iters);

	% Check that the provided w is valid
	if ~is_finite_acs(w(1:N) + 1j*[0;w(N+1:end)])
		error(['The provided w is not a valid finite autocorrelation sequence (it is not in C^*). ', ...
				'The problem is unbounded below.']);
	end

	% Initial values
	u = u0;

	% Initialize derivobj
	dobj = derivobj();
	if method == 2
		dobj.allow_form_hessian = true;
	end
	dobj.y = y;
	dobj.tau = tau;
	dobj.w = w;
	dobj = dobj.set_u(u);

	% Increase u(1) until we have found a primal-dual feasible starting point
	k = 0;
	while(true)
		k = k + 1;

		% Calculate duality gap to initilize t
		x = dobj.x;
		v = real(dobj.phi' * x);
		rho = tau;
		s = 2*(x-y);
		z = tau * w;
		% Check for dual feasibility
		in_dual = in_dual_cone(rho, s, z);
		if in_dual
			break;
		else
			u(1) = u(1) * 2;
			dobj = dobj.set_u(u);
		end
		if k > 100
			error('Could not find dual-feasible starting point');
		end
	end
	maxdual = -norm(s)^2/4-real(s'*y);
	gap = rho*v + z'*u + real(s'*x);
	dobj.t = mu*(N+1)/abs(gap);

	% Initialize L-BFGS object
	if method == 1
		lbfgs_obj = lbfgs_init(num_saved_lbfgs, 2*N-1);
	else
		lbfgs_obj = lbfgs_init(1, 1);	% Declared for codegen, unused
	end

	% Outer loop (varying barrier parameter)
	k = 0;
	gap = nan;
	grad1 = dobj.grad1;
	grad2 = dobj.grad2;
	while( true )
		k = k + 1;

		%%%%
		%% Obtain search direction
		%%%%
		grad = grad1 + 1/dobj.t * grad2;
		extra_info = '';
		if method == 1
			tmp = (N-1:-1:1)' / (2*N);
			diag_H = [1; tmp; tmp] * dobj.u0_hess;
			du = lbfgs_get_direction(lbfgs_obj, diag_H, grad, dobj.t);
			% If gradient is almost orthogonal to search direction, restart
			% L-BFGS
			if abs(du'*grad) / (norm(grad)*norm(du)) < 1e-5
				extra_info = ', restart';
				du = -grad;
				lbfgs_obj = lbfgs_restart(lbfgs_obj);
			end
		else
			du = - dobj.form_hessian \ grad;	% Exact Newton step
		end
		assert(~any(isnan(du)), 'du is nan.');
		if du'*grad >= 0
			if verbose
				fprintf('The search direction is not a descent direction. Using gradient\n');
			end
			du = -grad;
			if method == 1
				lbfgs_obj = lbfgs_restart(lbfgs_obj);
			end
		end


		%%%%
		%% Line search
		%%%%
		u0 = u;
		g0 = dobj.g;
		grad0 = grad;
		grad1_old = grad1;
		grad2_old = grad2;
		lsiters = 0;
		stepsize = 1;
		while true
			lsiters = lsiters + 1;
			% Find new values
			u = u0 + stepsize*du;
			[dobj, feas] = dobj.set_u(u);
			if all(abs(u-u0)<1e-10)
				% Progress is slow, force an increase of t
				[dobj, feas] = dobj.set_u(u0);
				dobj.t = mu*dobj.t;
				extra_info = [extra_info, ', inc t'];
				% Do not allow this to happen too many times
				if dobj.t > mu^3 * (N+1)/gap
					if verbose | print_warnings
						fprintf('Progress stagnated before acheiving requested tolerance.\n');
						fprintf('Terminating with current solution.\n');
					end
					x = dobj.x;
					v = 1/(tau*dobj.t) + real(dobj.phi' * x);
					return;
				else
					break;
				end
			end
			if feas
				grad1 = dobj.grad1;
				grad2 = dobj.grad2;
				grad = grad1 + 1/dobj.t * grad2;
				if dobj.g < g0 + 0.05*stepsize*(du'*grad0)
					break
				end
			end
			%if stepsize<1e-25
				%fprintf('Error: Line search did not find an acceptable stepsize. Exiting...\n');
				%lbfgs_obj = lbfgs_restart(lbfgs_obj);
				%[dobj, feas] = dobj.set_u(u0);
				%grad1 = grad1_old;
				%grad2 = grad2_old;
				%grad = grad0;
				%x = dobj.x;
				%return;
			%end
			
			stepsize = stepsize_reduction * stepsize;
		end

		%%%%
		%% Add new point to L-BFGS approximation
		%%%%
		if method==1
			yy1 = grad1 - grad1_old;
			yy2 = grad2 - grad2_old;
			ss = u - u0;
			lbfgs_obj = lbfgs_add_point(lbfgs_obj, ss, yy1, yy2);
		end

		%%%%
		%% Calculate duality gap
		%%%%
		% Get primal- and dual variables
		x = dobj.x;
		v = 1/(tau*dobj.t) + real(dobj.phi' * x);
		rho = tau;
		s = 2*(x-y);
		z = tau * w;
		objective = norm(x-y)^2 + tau*v + tau*w'*u;
		% Check for dual feasibility
		in_dual = in_dual_cone(rho, s, z);

		if in_dual
			if -norm(s)^2/4-real(s'*y) > maxdual
				saves = s;
			end
			maxdual = max(maxdual, -norm(s)^2/4-real(s'*y));
		end
		gap = objective - maxdual;
		relgap = gap / objective;
		if gap < 0
			if verbose | print_warnings
				fprintf('Encountered negative duality gap before acheiving requested tolerance.\n');
				fprintf('This is due to issues with numerical precision in assessing dual feasibility.\n');
				fprintf('Terminating with current solution.\n');
			end
			%if in_dual
				%assert( v*rho + real(x'*s) + z'*u < 0 );
			%end
			x = dobj.x;
			v = 1/(tau*dobj.t) + real(dobj.phi' * x);
			return;
		end

		% Save to info
		if isnan(gap) && k > 1
			info.gap(k) = info.gap(k-1);
			info.relgap(k) = info.relgap(k-1);
		else
			info.gap(k) = gap;
			info.relgap(k) = relgap;
		end
		info.objective(k) = objective;
		info.t(k) = dobj.t;
		info.g(k) = dobj.g;
		info.v(k) = v;
		info.x(:,k) = x;
		info.u(:,k) = u;
		info.iters = k;

		%%%%
		%% Print some info
		%%%%
		if verbose
			fprintf('%3i. gap=%7.1e / %7.1e  dual_feas=%i  t=%7.1e  grad=%7.1e  iters=(%i,%i)%s\n', ...
					int32(k), gap, relgap, int32(in_dual), dobj.t, ...
					norm(grad), int32(lsiters), ...
					int32(dobj.num_calls), extra_info);
		end


		% Stopping criterion
		if ~isnan(gap) && gap>=0
			if gap<tol_abs || relgap<tol_rel
				break;
			end
		end
		if k >= max_iters
			if verbose | print_warnings
				fprintf('Did not converge before reaching maximum number of iters.\n');
			end
			break;
		end

		% Update t
		if ~isnan(gap)
			dobj.t = max( dobj.t, mu*(N+1)/gap );
		end
	end

	% Get final x and v
	x = dobj.x;
	v = 1/(tau*dobj.t) + real(dobj.phi' * x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check for dual feasibility
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = in_dual_cone(rho, s, z)
	N = length(s);

	%vec = z - toeplitz_r2c_adjoint(s*s')/(4*rho);
	%vec = vec(1:N) + [0;1j*vec(N+1:end)];
	z_complex = z(1:N) + [0;1j*z(N+1:end)];
	FFTs = fft(s, 2*N-1);
	adj = 2*conj(ifft(abs(FFTs).^2));
	vec = z_complex - adj(1:N) / (4*rho);
	res = is_finite_acs(vec);

	if 0
		cvx_quiet true
		cvx_solver SDPT3
		cvx_begin
			variable T(N, N) hermitian toeplitz
			variable v(1)
			variable x(N) complex

			minimize( rho*v + real(s'*x) + z(1)*T(1,1)/2 )
			subject to
				[T x; x' v] == hermitian_semidefinite(N+1)
		cvx_end
		f = rho*v + real(s'*x) + z(1)*T(1,1)/2;
		fprintf('    CVX: %s, f=%.4g,  crit=%.4g\n', cvx_status, f, crit);
		if ~isempty(strfind(cvx_status, 'Solved'))
			res = true;
		else
			res = false;
		end
	end
end

% Is the complex vector z a valid finite autocorrelation sequence?
function res = is_finite_acs(z)
	N = length(z);
	Nfft = 2^ceil(log2(128*N));
	z(1) = z(1)/2;
	res = all( real(fft(z,Nfft))/sqrt(Nfft) >= 0 );
end

