% Copyright (c) 2018: Thomas Lundgaard Hansen.
% This software is being released under the MIT license (see LICENSE file). 
% Parse varagin
function as = parse_varargin(as_default, arglist)
	as = as_default;

	% No extra arguments given, use defaults
	if numel(arglist) == 0
		return;
	end

	% Handle both structs and <name>,<value> pairs
	if numel(arglist)==1 & isstruct(arglist{1})
		input = arglist{1};
	else
		input = struct(arglist{:});
	end

	% Loop over names in input and assign given values to parm
	names = fieldnames(input);
	for i = 1:length(names)
		name = names{i};
		if isfield(as, name)
			as.(name) = input.(name);
		else
			error('No field in as_default named ''%s''', name);
		end
	end
end
