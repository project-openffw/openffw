function p = computeSolution(p)
% Main computation loop. 
% Computes the AFEM loop until stopping criterion is
% reached and specifies functions to evaluate the error

% Copyright 2007 Jan Reininghaus, David Guenther, 
%                Andreas Byfut, Joscha Gedicke
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


%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set initial level number
curLevel = 1;
p.level(1).level = curLevel;

% load abort criteria
maxLevel = loadField('p.params','maxLevel',p,1000);
minLevel = loadField('p.params','minLevel',p,1);
maxNrDoF = p.params.maxNrDoF;

% refine first levels
refineFirstLevels = loadField('p.params.modules.mark','refineFirstLevels',p,1);

% handles for first level calculations
mark = @uniform;
refine = p.statics.refine;
solve = p.statics.solve;
createLinSys = p.statics.createLinSys;
estimate = p.statics.estimate;
enumerate = p.statics.enumerate;
postProc = p.statics.postProc;

%% Initial enumeration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = enumerate(p);

%% Refine First Level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to refine the initial mesh uniformly as often as
% spezified in 'refineFirstLevels'.

if refineFirstLevels > 0
    q = p;
    while( curLevel <= refineFirstLevels )
        q = mark(q);
        q = refine(q);
        q = enumerate(q);
        curLevel = curLevel+1;
        q.level(curLevel).level = curLevel;
    end
    curLevel = 1;
    p.level = q.level(end);
    p.level(1).level = curLevel;
    clear q;  % clear memory
end

%% Compute Solution for the First Level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% builds matrix A and RHS b
p = createLinSys(p);

% computes the solution of A*x = b
p = solve(p);

% post processes the solution x
p = postProc(p); 

% estimate error
p = estimate(p);


%% AFEM loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n');

nrDoF = p.level(end).nrDoF;
fprintf(stderr,'Level = %d \t Error = %.2g \t DoF = %d\n',curLevel,p.level(end).estimatedError,nrDoF);
while(nrDoF < maxNrDoF && curLevel < maxLevel) || (curLevel < minLevel)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = afem(p);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    curLevel = curLevel + 1;
    nrDoF = p.level(end).nrDoF;
    estimatedError = p.level(end).estimatedError;
    fprintf(stderr,'Level = %d \t Error = %.2g \t DoF = %d\n',curLevel,estimatedError,nrDoF);
    p.level(end).level = curLevel;    
    maxNrDoF = p.params.maxNrDoF;

end



