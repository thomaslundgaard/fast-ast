function out = trTinvTinv(rho, delta)
	%#codegen
	N = length(rho);
	i_sum = zeros(N,1);
	rho_zero = [0; rho];

	for i = 1:N
		mn_sum = conj(rho(N:-1:i)).*rho(N-(i-1):-1:1) ...
			- rho_zero(1:N-i+1).*conj(rho_zero(i:N));
		mn_sum = cumsum(mn_sum);
		mn_sum = abs(mn_sum).^2;
		i_sum(i) = 2 * sum(mn_sum);
	end
	i_sum(1) = i_sum(1) / 2;	% undo scaling by 2
	out = sum(i_sum);
	out = out / delta^2;
end

