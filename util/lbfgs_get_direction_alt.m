% Copyright (c) 2018: Thomas Lundgaard Hansen.
% This software is being released under the MIT license (see LICENSE file). 
function delta = lbfgs_get_direction(obj, initial_diag_H, deriv, t)
	%#codegen
	num_saved = obj.num_saved;
	assert(num_saved <= size(obj.ss,2));	% Allow static memory allocation

	H = diag( initial_diag_H );

	% Find index of oldest saved vector
	k = mod( (obj.newest+num_saved-1)-1, obj.max_saved ) + 1;
	for i = 1:num_saved
		yy = obj.yy1(:,k) + 1/t * obj.yy2(:,k);
		Hs = H * obj.ss(:,k);
		H = H + yy*yy' / (yy'*obj.ss(:,k)) - Hs * Hs' / (obj.ss(:,k)'*Hs);
		k = mod( (k-1)-1, obj.max_saved ) + 1;
	end

	% Find direction
	delta = - H \ deriv;
end

