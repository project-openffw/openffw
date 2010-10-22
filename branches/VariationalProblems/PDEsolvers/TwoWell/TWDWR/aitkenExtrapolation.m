function p = aitkenExtrapolation(value,p)

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


s = zeros(length(p.level),1);
for curLvl = 2 : length(p.level)
    s(curLvl) = getfield(p.level, {curLvl} , value);
end

while length(s) >= 3
    N = length(s);
    s = s(1:N-2)- (s(2:N-1)-s(1:N-2)).^2./(s(3:N)-2*s(2:N-1)+s(1:N-2));
end

p = setfield(p, ['aitkenExtrapolation',value], s(end));

% for curLvl = 1 : length(p.level)
%     t = getfield(p.level, {curLvl} , value);
%     p = setfield(p, ['level(',int2str(curLvl),').',value,'Error'], norm(t-s(end)));
% end

