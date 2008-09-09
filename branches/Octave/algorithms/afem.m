function p = afem(p) 
% AFEM loop calculations
% i.e., mark => refine => enumerate => createLinSys => solve 
%       => postProc => estimate

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


%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mark = p.statics.mark;
refine = p.statics.refine;
solve = p.statics.solve;
createLinSys = p.statics.createLinSys;
estimate = p.statics.estimate;
enumerate = p.statics.enumerate;
postProc = p.statics.postProc;


%% AFEM loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mark edges for refinement
p = mark(p);

% RedGreenBlue refinement of markededges
p = refine(p);

% create necessary data structures
p = enumerate(p);

% builds matrix A and RHS b
p = createLinSys(p);

% computes the solution of A*x = b
p = solve(p);

% post processes the solution x
p = postProc(p); 

% estimate error
p = estimate(p);
