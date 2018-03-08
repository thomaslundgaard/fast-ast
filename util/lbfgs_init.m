% Copyright (c) 2018: Thomas Lundgaard Hansen.
% This software is being released under the MIT license (see LICENSE file). 
function obj = lbfgs_init(max_saved, N)
	obj = struct();
	obj.max_saved = int32(max_saved);
	obj.num_saved = int32(0);
	obj.newest = int32(1);
	obj.ss = nan(N, max_saved);
	obj.yy1 = nan(N, max_saved);
	obj.yy2 = nan(N, max_saved);
	obj.sy1 = nan(max_saved, 1);
	obj.sy2 = nan(max_saved, 1);
end

