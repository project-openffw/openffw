function p = graded(p)
% graded meshes
% generates graded meshes using the red-green-blue refinement

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

level = p.level(end).level;
beta = loadField('p.params.modules.mark.graded','beta',p,1/3);
N = loadField('p.params.modules.mark.graded','gradeN',p,1);


%% Graded Refinement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = (1/2)^N;
% max_h = (6/7)^N;
mu = beta;

count = 1;
nextLevel = level + 1;

p.level(nextLevel).geom = p.level(nextLevel-1).geom;
enumerate = p.statics.enumerate;
refine = p.statics.refine;
 p = enumerate(p);
% 	p.level(nextLevel-1).enum.newNode4ed = [];

while count ~= 0
    count = 0;
    
    ed4e = p.level(end).enum.ed4e;
    midPoint4e = p.level(end).enum.midPoint4e;
    area4e = p.level(end).enum.area4e;
    nrEdges = p.level(end).nrEdges;
    nrElems = p.level(end).nrElems;
    
    markedEdges = false(nrEdges,1);
    refineElems = false(nrElems,1);
    
    for curElem = 1:nrElems
        mP = midPoint4e(curElem,:);
        area = area4e(curElem,:);
        hT = sqrt(area);
        PhiMuT = norm(([0,0] - mP))^(1-mu);
        if hT > H*PhiMuT;
            refineElems(curElem) = true;
            count = count + 1;
        end
    end

    markedEdges4e = ed4e(refineElems,:);
    markedEdges(markedEdges4e(:)) = true;
    
    p.level(end).markedEdges = markedEdges';
    q = p.statics.closure(p);
    q.level(end).refineElemsBisec5 = false(nrElems,1);
    q.level(end).level = length(q.level);
    q = refine(q);
    p.level(end).geom = q.level(end).geom;
    p = enumerate(p);
end 

nrEdges = p.level(end).nrEdges;
nrElems = p.level(end).nrElems;
markedEdges = false(nrEdges,1);

%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.level(end).markedEdges = markedEdges';
p.level(end).refineElemsBisec5 = false(nrElems,1);
p = p.statics.closure(p);
p.params.modules.mark.graded.gradeN = N+1;
