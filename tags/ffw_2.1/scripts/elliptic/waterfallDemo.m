function p = waterfallDemo
% waterfall example

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


%% Problem Definition
problem = 'Elliptic_Waterfall_exact';

pdeSolver = 'P1';         %'P1-Elliptic';
% pdeSolver = 'CR';         %'CR-Elliptic';
% pdeSolver = 'RT0P0';      %'RT0-P0-mixed-Elliptic';

maxNrDoF = 100;

% drawFunc = 'drawError_estimatedError';
drawFunc = 'drawError_H1semiError';
% drawFunc = 'drawError_L2error';

%% Patch Run
K = 1:10:51;
level = 1;
quotient = zeros(1,length(K));

for curK = K

	mark = 'bulk';
    p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elliptic');
    p.PDE.K = curK;
    p = computeSolution(p);
    
    p.params.output.holdIt = true;
    p.params.output.name = [mark,', ','K = ',int2str(curK)];
    p = show(drawFunc,p);

	nrDoF = p.level(end).nrDoF;
	bulkError = p.level(end).estimatedError;
	bulkError = bulkError/sqrt(nrDoF);

	mark = 'uniform';
    p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elliptic');
    p.PDE.K = curK;
    p = computeSolution(p);

    p.params.output.holdIt = true; 
    p.params.output.name = [mark,', ','K = ',int2str(curK)];
    p = show(drawFunc,p);

	nrDoF = p.level(end).nrDoF;
	unifError = p.level(end).estimatedError;
	unifError = unifError/sqrt(nrDoF);

	quotient(level) = unifError/bulkError;
	level = level +1;
end

figure(2)
plot(K,quotient)

