function regTestEllipticPDEs
% test script for all implemented elliptic problems

% Copyright 2007 Jan Reininghaus
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


maxNrDoF = 100;


pdeSolvers = struct('z',{'RT0P0','P1','CR'} );				
problems = struct('z',{'Elliptic_Square_exact',...
					   'Elliptic_HexagonalSlit_exact',...
                       'Elliptic_Lshape_exact',...
                       'Elliptic_Waterfall_exact'} );
                  
errTypes = struct('z',{'estimatedError','H1semiError','L2error'} );
marking = struct('z',{'bulk','uniform'});
				  
for i = 1:length(pdeSolvers)
    pdeSolver = pdeSolvers(i).z;
    for j = 1:length(problems)
        figure
        problem = problems(j).z;
        for m = 1:length(marking)
            mark = marking(m).z;
            p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elliptic');
            p = computeSolution(p);
            for k = 1:length(errTypes)
                errType = errTypes(k).z;

                drawErrorString = ['drawError_',errType];
                p.params.output.drawConvergenceRate = true;
                p.params.errorIntegtrateExactDegree = 9;
                p.params.output.holdIt = true;
                p.params.output.drawGrid = false;
                p = show(drawErrorString,p);

            end
        end
    end
end
		
	
