% Copyright (c) 2018: Thomas Lundgaard Hansen.
% This software is being released under the MIT license (see LICENSE file). 
function obj = lbfgs_restart(obj)
	obj.num_saved = int32(0);
	obj.newest = int32(1);
end

