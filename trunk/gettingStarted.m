function p = gettingStarted
% First start with the FFW

% Copyright 2007 Andreas Byfut
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


problem = 'Elliptic_Lshape';
pdeSolver = 'P1';
maxNrDoF = 100;
mark = 'bulk';

p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elliptic');
p = computeSolution(p);

figure(1);
set(gcf,'Name','Displacement');
p = show('drawU',p);
view(-30,20)

figure(2);
set(gcf,'Name','Estimated Error');
p = show('drawError_estimatedError',p);
