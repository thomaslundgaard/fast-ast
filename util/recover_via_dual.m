function omega = recover_via_dual(x, y, tau)
	% Inspired by dual_poly_debias.m from
	% https://github.com/badrinarayan/astlinespec

	N = length(y);

	grid_size = 2^ceil(log2(2048*N));

	amp = grid_size / tau; % adjust dual polynomial amplitude to peak at 1

	dual_poly_coeff = (y - x);

	polyval = amp * abs(ifft(conj(dual_poly_coeff(:)), grid_size));

	%hold off; plot( linspace(0,1-1/grid_size,grid_size), polyval);

	polyval(polyval < 1-2e-2) = 0;

	%hold on; plot(linspace(0,1-1/grid_size,grid_size), polyval);

	[~, omega] = findpeaks(polyval);

	if polyval(1)>polyval(2) && polyval(1)>polyval(end)
		omega = [1;omega];
	end
	if polyval(end)>polyval(1) && polyval(end)>polyval(end-1)
		omega = [omega; grid_size];
	end

	omega = (omega(:)-1) / grid_size;
end

