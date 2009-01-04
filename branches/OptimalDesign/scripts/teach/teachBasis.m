function teachBasis
% shows the shape of the basis functions

% Copyright 2007 Jan Reininghaus, David Guenther
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


problem = 'Elliptic_Square';
mark = 'uniform';
maxNrDoF = 1;

%% PDE Solver
% pdeSolver = 'P1';
% pdeSolver = 'P2';
pdeSolver = 'P3';
% pdeSolver = 'CR';

%% Init
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elliptic');
p = p.statics.enumerate(p);

%% Draw Basis
nrBasisFunc = size(p.level(end).enum.dofU4e,2);
p.params.output.holdIt = false;

for curBasisFunc = 1 : nrBasisFunc
    drawDisplacementBasis(curBasisFunc,100,p);
    view(5,20)
    shading interp
    lighting gouraud
    light('Position',[0 0 2],'Style','infinite');
    pause
end
