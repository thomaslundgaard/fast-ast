% Copyright (c) 2018: Thomas Lundgaard Hansen.
% This software is being released under the MIT license (see LICENSE file). 
function obj = lbfgs_add_point(obj, ss, yy1, yy2)
	obj.newest = mod( (obj.newest-1)-1, obj.max_saved ) + 1;
	obj.num_saved = min(obj.num_saved+1, obj.max_saved);
	ss(ss==0) = eps();
	yy1(yy1==0) = eps();
	yy2(yy2==0) = eps();
	obj.ss(:, obj.newest) = ss;
	obj.yy1(:, obj.newest) = yy1;
	obj.yy2(:, obj.newest) = yy2;
	obj.sy1(obj.newest) = ss' * yy1;
	obj.sy2(obj.newest) = ss' * yy2;
end

