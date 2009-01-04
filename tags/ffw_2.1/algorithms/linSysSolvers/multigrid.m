function p = multigrid(p)
% Solves the system 'Ax=b' iteratively.
% The iterative multigrid solver, using Richardson's smoother.

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
A = p.level(end).A;
b = p.level(end).b;
x = loadField('p.level(end)','x',p);

fixedNodes = p.level(end).enum.fixedNodes;
freeNodes = p.level(end).enum.freeNodes;
markedEdges = loadField('p.level(end-1)','markedEdges',p);

damp = loadField('p.params.modules.solver.multigrid','damp',p,0.2);
nested = loadField('p.params.modules.solver.multigrid','nested',p,true);
maxSteps = loadField('p.params.modules.solver.multigrid','maxSteps',p,1000);


%% MULTIGRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AFree = A(freeNodes,freeNodes);
bFree = b(freeNodes);

% Solve the system on the coarsest grid.
if( p.level(end).level == 1)
	x(freeNodes) = AFree \ bFree;
	prolong = eye(length(freeNodes));
	p.level(1).solve.prolong = prolong;
    
% Start of multigrid cycle.
else
    % get information of former grids
    xOld = loadField('p.level(end-1)','x',p);
    fixedNodesOld = p.level(end-1).enum.fixedNodes;
    markededges = loadField('p.level(end-1)','refineEdges',p);
    nrNodes = p.level(end).nrNodes;
    nrNodesOld = p.level(end-1).nrNodes;
    n4edOld = loadField('p.level(end-1).enum','n4ed',p);
    newNode4ed = loadField('p.level(end-1).enum','newNode4ed',p);
    estError = p.level(end-1).estimatedError;
    
   	oldDim = nrNodesOld;
	dim = nrNodes;

    % Create prolongation matrix
	prolong = spdiags( ones(oldDim,1),0,oldDim,oldDim);
	dummy = [newNode4ed(markedEdges); newNode4ed(markedEdges)];
	I = dummy(:) - oldDim;
	dummy = n4edOld(markedEdges,:)';
	J = dummy(:);
	S = 1/2 * ones(2*(dim - oldDim),1);
	weights = sparse(I,J,S, (dim - oldDim), oldDim);

	prolong = [prolong; weights];

	xOld = prolong*xOld;

	prolong(fixedNodes,:) = [];
	prolong(:,fixedNodesOld) = [];
    
    % Save prolongation matrix
	p.level(end-1).solve.prolong = prolong;

	lvl = p.level(end-1).level + 1;
	if(nested)
		xFree = xOld(freeNodes);
	else
		xFree = zeros(length(freeNodes),1);
	end
% 	tol = (1e-4)/sqrt(length(freeNodes));
    tol = estError/1e3;
	eta = 1e10;

    % Start the multigrid recursion
	steps = 0;
	while eta > tol
		resid = bFree - AFree*xFree;
		etaOld = sqrt(resid'*AFree*resid);
		steps = steps + 1;
        
        % Start Multigrid
		xFree = mg(xFree,lvl,bFree,p,damp);
        
		resid = bFree - AFree*xFree;
		eta = sqrt(resid'*AFree*resid);
		if(eta/etaOld > 1)
			damp = damp / 1.5;
		end
		if(steps > maxSteps)
			error('unstable damping parameter');
        end
  		% rate = eta/etaOld;
    end
	x(freeNodes) = xFree;
end


%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.level(end).x = x;


%% Internal Multigrid Steps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = mg(x0,lvl,b,p,damp)
% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prolong = p.level(lvl-1).solve.prolong;
A = p.level(lvl).A;
freeNodes = p.level(lvl).enum.freeNodes;
A = A(freeNodes,freeNodes);
preSmoothSteps = loadField('p.params.modules.solver.multigrid','preSmoothSteps',p,3);
postSmoothSteps = loadField('p.params.modules.solver.multigrid','postSmoothSteps',p,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(lvl == 2)
    % direct solve if we are on lowest level
	x = A \ b;
else
    % pre Smoothing
	x = richard(A,b,x0,damp,preSmoothSteps);

	restrict = prolong';

	defect = b - A*x;
	defect = restrict*defect;

	x0 = zeros(length(defect),1);
	correction = mg(x0,(lvl-1),defect,p,damp);

	correction = prolong*correction;
	x = x + correction;

    % post Smoothing
    x = richard(A,b,x,damp,postSmoothSteps);
end


%% The Richardson Smoother %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = richard(A,b,x,damp,nrSteps)
for curStep = 1:nrSteps
	x = x - damp * (A*x-b);
end
