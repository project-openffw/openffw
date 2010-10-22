function [a,b] = binsearch(X,v)

% Copyright 2007 Joscha Gedicke
%
% This file is part of FFW.
%
% FFW is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% FFW is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

N = length(X);

index = floor(N/2);
a=1;
b=N;
w = X(index);

while( w ~= v )
    if w > v
        b=index;
    else
        a=index;
    end
    index = a+ceil((b-a)/2);
    w = X(index);
end

a=index;
while (a>1) && X(a-1)==v
    a = a-1;
end

b=index;
while (b<N) && X(b+1)==v
    b = b+1;
end
