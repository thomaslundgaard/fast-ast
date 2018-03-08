% Copyright (c) 2018: Thomas Lundgaard Hansen, Tobias Lindstrøm Jensen.
% This software is being released under the MIT license (see LICENSE file). 
classdef derivobj
% Handles calculation of our objective, gradients thereof and the primal
% variables.

properties
	y
	tau
	w
	t
	num_calls = 0;
	allow_form_hessian = false;
end % properties

properties (Constant)
	check_against_direct = false;
end % properties

properties (SetAccess = private)
	u = nan;
	FFTy = nan;
	rho = nan;
	R = nan;
	delta = nan;
	delta_tau = nan;
	FFTrho_tau = nan;
	FFTflipconjrho_tau = nan;
	phi = nan;
	FFTphi = nan;
	% Used by direct methods
	Tu
	TuInv
	Tu_tauInv
end % properties

properties (Dependent = true)
	grad1	% real-valued
	grad2	% real-valued
	x		% complex-valued
	g
end % properties

methods
	%%%%%%%%%%%%%%%%%%%%
	%% Set y and u
	%%%%%%%%%%%%%%%%%%%%
	function obj = set.y(obj, y)
		N = length(y);
		% Find FFT size as smallest power-of-two which doesnt give aliasing
		Nfft = 2^ceil(log2(2*N - 1));
		obj.FFTy = fft(y, Nfft);
		obj.y = y;
	end

	function [obj, feasible] = set_u(obj, u)
		assert(iscolumn(u));
		N = (length(u)+1) / 2;
		Nfft = length(obj.FFTy);

		% Fast Toeplitz inversion of T(u)
		row = [2*u(1); u(2:N) + 1j*u(N+1:end)];	% First row of T(u)
		row = row(:).';	% Make sure its a row vector
		if obj.allow_form_hessian
			[obj, rho, delta, R] = obj.fast_toeplitz_inv(row);
		else
			[obj, rho, delta] = obj.fast_toeplitz_inv(row);
			R = nan;	% For codegen
		end

		% Check feasibility
		feasible = all(delta > 0);
		if obj.check_against_direct
			feasible_old = all(eig(toeplitz_r2c(u))>100*eps());
			assert(feasible == feasible_old);
		end
		if ~feasible 
			if nargout < 2
				fprintf('Warning: Infeasible u and ''feasible'' output is ignored.\n');
			end
			return;
		end

		% Fast Toeplitz inversion of T(u + tau/2 e_0)
		row(1) = row(1) + obj.tau;
		[obj, rho_tau, delta_tau] = obj.fast_toeplitz_inv(row);

		% Calculate phi
		FFTrho_tau = fft(rho_tau, Nfft);
		FFTflipconjrho_tau = fft(conj(flipud(rho_tau)), Nfft);
		phi = obj.toeplitz_inv_prod(FFTrho_tau, FFTflipconjrho_tau, delta_tau, obj.FFTy);

		% Assign values
		obj.u = u;
		obj.rho = rho;
		if obj.allow_form_hessian
			obj.R = R;
		end
		obj.delta = delta;
		obj.delta_tau = delta_tau;
		obj.FFTrho_tau = FFTrho_tau;
		obj.FFTflipconjrho_tau = FFTflipconjrho_tau;
		obj.phi = phi;
		obj.FFTphi = fft(phi, Nfft);

		% Direct calculation
		if obj.check_against_direct
			obj.Tu = toeplitz_r2c(obj.u);
			obj.TuInv = inv(obj.Tu);
			u_tau = u;
			u_tau(1) = u_tau(1) + obj.tau/2;
			obj.Tu_tauInv = inv(toeplitz_r2c(u_tau));
			phi_old = obj.Tu_tauInv * obj.y;
			obj.assert_equal(obj.phi, phi_old);
		end
	end

	%%%%%%%%%%%%%%%%%%%%
	%% Calculate objective and gradients
	%%%%%%%%%%%%%%%%%%%%
	function out = get.g(obj)
		out = obj.tau*obj.w'*obj.u ...
				+ obj.tau * real(obj.y' * obj.phi) ...
				- 1/obj.t * sum(log(obj.delta));

			if obj.check_against_direct
			old = obj.tau*obj.w'*obj.u ...
				+ obj.tau * real(obj.y' * obj.phi) ...
				- 1/obj.t * logdet(obj.Tu, 'chol');
			obj.assert_equal(out, old);
		end
	end

	function grad1 = get.grad1(obj)
		N = length(obj.y);
		tmp = - 2 * conj(ifft(abs(obj.FFTphi).^2));
		tmp = tmp(1:N);
		tmp = [real(tmp); imag(tmp(2:end))];

		if obj.check_against_direct
			tmp_old = - toeplitz_r2c_adjoint(obj.phi*obj.phi');
			obj.assert_equal(tmp, tmp_old);
		end

		grad1 = obj.tau * obj.w + obj.tau * tmp;
	end

	function grad2 = get.grad2(obj)
		% Calc omega_s based on method in T.L. Hansen et al - "Superfast Line
		% Spectral Estimation", IEEE Trans. Signal Process.
		N = length(obj.y);
		Nfft = length(obj.FFTy);
		FFTrho = fft(obj.rho, Nfft);
		qq = [0:N-1]';
		omega1 = conj(ifft( abs(FFTrho).^2 ));
		omega1 = (2-N+qq) .* omega1(1:N);
		tmp = fft(qq.*obj.rho, Nfft);
		omega2 = conj(ifft( conj(tmp) .* FFTrho ));
		omega2 = 2*omega2(1:N);
		omega_s = (omega1 + omega2) / obj.delta(end);

		grad2 = - 2 * omega_s;
		grad2 = [real(grad2); imag(grad2(2:end))];

		if obj.check_against_direct
			grad2_old = - toeplitz_r2c_adjoint(obj.TuInv);
			obj.assert_equal(grad2, grad2_old);
		end
	end

	%%%%%%%%%%%%%%%%%%%%
	%% Calculate x
	%%%%%%%%%%%%%%%%%%%%
	function x = get.x(obj)
		N = length(obj.y);

		% Get the conjugate-symmetric complex version of u (i.e. the autocorrelation)
		tmp = obj.u(2:N)+1j*obj.u(N+1:end);
		uu = [conj(flipud(tmp)); 2*obj.u(1); tmp];
		tmp = ifft(fft(conj(uu),length(obj.FFTphi)) .* obj.FFTphi);
		x = tmp(N:2*N-1);

		if obj.check_against_direct
			x_old = obj.Tu * obj.phi;
			obj.assert_equal(x, x_old);
		end
	end

	%%%%%%%%%%%%%%%%%%%%
	%% Calculate values related to the Hessian
	%%%%%%%%%%%%%%%%%%%%
	function out = u0_hess(obj)
		N = length(obj.y);
		Nfft = length(obj.FFTy);

		% Calculate first term: phi * Tu_tauInv * phi
		T0hphi = ifft( obj.FFTflipconjrho_tau .* obj.FFTphi );
		T0hphi = T0hphi(N+1:2*N);
		T0hphi(end) = 0;
		T1phi = ifft( obj.FFTrho_tau .* obj.FFTphi );
		T1phi = T1phi(N:2*N-1);
		phi_Tu_tauInv_phi = (norm(T1phi)^2 - norm(T0hphi)^2 ) ...
								/ obj.delta_tau(end);
		
		out1 = 8 * obj.tau * phi_Tu_tauInv_phi;

		% Calculate second term
		if isreal(obj.rho)
			rho = complex(obj.rho);
		else
			rho = obj.rho;
		end
		out2 = 4 / obj.t * trTinvTinv(rho, obj.delta(end));
		out = out1 + out2;

		if obj.check_against_direct
			out1 = 8 * obj.tau * real(obj.phi' * obj.Tu_tauInv * obj.phi);
			out2 = 4 / obj.t * sum(sum(abs(obj.TuInv).^2));
			out_old = out1 + out2;
			obj.assert_equal(out, out_old);
		end
	end

	function H = form_hessian(obj)
		N = length(obj.y);
		Nfft = length(obj.FFTy);
		H = zeros(2*N-1,2*N-1);
		if obj.check_against_direct
			H_barrier_old = zeros(2*N-1,2*N-1);
		end
		for i = 1:2*N-1
			%% Calculate a column of Hessian of g
			% Form a
			if i<=N
				shift = i-1;
				coef = 1;
			else
				shift = i-N;
				coef = 1j;
			end
			a = [zeros(shift,1); conj(coef)*obj.phi(1:N-shift)] ...
				+ coef*[obj.phi(shift+1:N); zeros(shift,1)];
			a = obj.toeplitz_inv_prod(obj.FFTrho_tau, obj.FFTflipconjrho_tau, ...
					obj.delta_tau, fft(a,Nfft));
			FFTa = fft(a, Nfft);

			% Cross-correlate with phi to get column of Hessian
			tmp = conj(ifft(2*real(obj.FFTphi.*conj(FFTa))));
			tmp = tmp(1:N);
			out = obj.tau * [2*real(tmp); 2*imag(tmp(2:end))];

			%% Check against direct calculation
			if obj.check_against_direct
				du = [zeros(i-1,1); 1; zeros(2*N-1-i, 1)];
				toeplitz_du = toeplitz_r2c(du);
				dphi = - obj.Tu_tauInv * (toeplitz_du*obj.phi);
				dphi_phi = dphi * obj.phi';
				out_old = - obj.tau * toeplitz_r2c_adjoint(dphi_phi+dphi_phi');
				obj.assert_equal(out, out_old);

				H_barrier_old(:,i) = toeplitz_r2c_adjoint(obj.TuInv * toeplitz_du * obj.TuInv);
			end

			H(:,i) = out;
		end

		% Form Hessian of the barrier function
		H_barrier = form_barrier_hessian(obj.R, obj.delta);
		if obj.check_against_direct
			obj.assert_equal(H_barrier, H_barrier_old);
		end

		% Calculate complete Hessian
		H = H + 1/obj.t * H_barrier;
	end

	%%%%%%%%%%%%%%%%%%%%
	%% Utility functions
	%%%%%%%%%%%%%%%%%%%%
	function [obj, rho, delta, R] = fast_toeplitz_inv(obj, row)
		obj.num_calls = obj.num_calls + 1;
		if length(row)<=512 || nargout>=4
			% Use levinson-durbin recursion
			if isreal(row)
				% Old versions are a bit stupid and need this spelled out
				r = complex(conj(row(:)));
			else
				r = conj(row(:));
			end
			if nargout>=4
				[a, P, R] = levin(r);
			else
				[a, P] = levin(r);
			end
			rho = flipud(a);
			delta = P;
		else
			% Use Generalized Schur algorithm
			[xi, eta, gamma] = genschur(-row(2:end), row(1:end-1));
			rho = [0; flipud(eta(:))] + [flipud(xi(:)); 0];
			delta = cumprod([real(row(1)); 1-abs(gamma(:)).^2]);	% Form delta
		end
	end

	function out = toeplitz_inv_prod(obj, FFTrho, FFTflipconjrho, delta, FFTvec)
		% Calculates the product of a vector onto a toeplitz inverse, i.e.,
		% inv(T(u)) * vec.
		% The Gohberg-Semencul factorization of inv(T(u)) must be given.

		Nfft = length(FFTvec);
		N = length(delta);

		% Notation: T0 = L; T1 = U
		T0h_vec = ifft( FFTflipconjrho .* FFTvec );
		T0h_vec = T0h_vec(N+1:2*N);
		T0h_vec(end) = 0;
		T0T0h_vec = [0; ifft( FFTrho .* fft(T0h_vec, Nfft) )];
		T0T0h_vec = T0T0h_vec(1:N);

		T1_vec = ifft( FFTrho .* FFTvec );
		T1_vec = T1_vec(N:2*N-1);
		T1hT1_vec = ifft( FFTflipconjrho .* fft(T1_vec, Nfft) );
		T1hT1_vec = T1hT1_vec(1:N);
		
		out = (T1hT1_vec - T0T0h_vec) / delta(end);
	end

	function assert_equal(obj, a, b)
		assert( max(abs(a(:)-b(:))) < 1e-7 * max(abs(b(:))) );
	end

end % methods
end % classdef

