function p = multigrid(p)

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = p.level(end).A;
b = p.level(end).b;
x = loadField('p.level(end)','x',p);
xOld = loadField('p.level(end-1)','x',p);

fixedNodes = p.level(end).enum.fixedNodes;
fixedNodesOld = p.level(end-1).enum.fixedNodes;
freeNodes = p.level(end).enum.freeNodes;
markededges = loadField('p.level(end-1)','refineEdges',p);

nrNodes = p.level(end).nrNodes;
nrNodesOld = p.level(end-1).nrNodes;
n4edOld = loadField('p.level(end-1).enum','n4ed',p);
newNode4ed = loadField('p.level(end-1).enum','newNode4ed',p);

estError = p.level(end-1).estimatedError;

damp = loadField('p.params.modules.solver.multigrid','damp',p,0.2);
nested = loadField('p.params.modules.solver.multigrid','nested',p,true);
maxSteps = loadField('p.params.modules.solver.multigrid','maxSteps',p,1000);


AFree = A(freeNodes,freeNodes);
bFree = b(freeNodes);

if( p.level(end-1).level == 1)
	x(freeNodes) = AFree \ bFree;
	prolong = eye(length(freeNodes));
	p.level(end-1).solve.prolong = prolong;
else
	oldDim = nrNodesOld;
	dim = nrNodes;

    % Create prolongation matrix
	prolong = spdiags( ones(oldDim,1),0,oldDim,oldDim);
	dummy = [newNode4ed(markededges); newNode4ed(markededges)];
	I = dummy(:) - oldDim;
	dummy = n4edOld(markededges,:)';
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
			damp = damp / 1.5
		end
		if(steps > maxSteps)
			error('unstable damping parameter');
		end
		rate = eta/etaOld;
    end
    eta;
    steps;
	x(freeNodes) = xFree;
end


%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).x = x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = mg(x0,lvl,b,p,damp)

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prolong = p.level(lvl-1).solve.prolong;
A = p.level(lvl).A;
freeNodes = p.level(lvl).enum.freeNodes;
A = A(freeNodes,freeNodes);
preSmoothSteps = loadField('p.params.modules.solver.multigrid','preSmoothSteps',p,3);
postSmoothSteps = loadField('p.params.modules.solver.multigrid','postSmoothSteps',p,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = richard(A,b,x,damp,nrSteps)

for curStep = 1:nrSteps
	x = x - damp * (A*x-b);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







