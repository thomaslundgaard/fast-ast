% Copyright (c) 2018: Thomas Lundgaard Hansen.
% This software is being released under the MIT license (see LICENSE file). 
function delta = lbfgs_get_direction(obj, initial_diag_H, deriv, t)
	%#codegen
	num_saved = obj.num_saved;
	assert(num_saved <= size(obj.ss,2));	% Allow static memory allocation
	alpha = nan(num_saved,1);
	rho = nan(num_saved,1);
	qq = deriv;
	k = obj.newest-1;	% index of current set of saved vectors
	for i = 1:num_saved
		k = mod( (k+1)-1, obj.max_saved ) + 1;
		rho(i) = 1 ./ (obj.sy1(k) + 1/t * obj.sy2(k));
		alpha(i) = rho(i) * (obj.ss(:,k)'*qq);
		qq = qq - alpha(i) * (obj.yy1(:,k) + 1/t*obj.yy2(:,k));
	end
	zz = qq ./ initial_diag_H;
	for i = num_saved:-1:1
		yy = obj.yy1(:,k) + 1/t * obj.yy2(:,k);
		beta = rho(i) * (yy'*zz);
		zz = zz + obj.ss(:,k) * (alpha(i)-beta);
		k = mod( (k-1)-1, obj.max_saved ) + 1;
	end
	delta = - zz;
end

