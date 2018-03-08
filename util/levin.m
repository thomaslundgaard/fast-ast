% The MIT License (MIT)
% 
% Copyright (c) 2015 Thomas Lundgaard Hansen
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.


% Levinson Recursion
% Gives forward linear prediction coefficients in 'a'
% and FLP mean square errors in P
% See Table 7.3 in [Manolakis - Statistical and Adaptive Signal Processing].
function [a, P, R] = levin(r)
%#codegen
	% Preallocate
	a = complex(zeros(length(r),1));
	P = zeros(length(r),1);
	if nargout >= 3
		R = complex(zeros(length(r),length(r)));
	end

	% Init
	P(1) = real(r(1));
	a(1) = 1;
	if nargout >= 3
		R(1,1) = a(1);
	end

	% Iterate
	for m = 2:length(r)
		if P(m-1) <= 0
			return;
		end
		beta = a(2:m-1).'*flipud(conj(r(2:m-1))) + conj(r(m));
		k = - beta/P(m-1);
		a(2:m-1) = a(2:m-1) + flipud(conj(a(2:m-1)))*k;
		a(m) = k;

		P(m) = P(m-1) + real(beta * conj(k));
		if nargout >=3
			R(1:m,m) = flipud(a(1:m));
		end
	end
end

