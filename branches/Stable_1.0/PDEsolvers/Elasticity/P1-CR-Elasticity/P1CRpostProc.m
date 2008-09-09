function p = P1CRpostProc(p)
%postproc.m computes the post-processing datas 
%for the Kouhia-Stenberg FE in linear elasticity. 
%
%authors: David Guenther, Jan Reininghaus
%%%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%%
n4e = p.level(end).geom.n4e;

lambda = p.PDE.lambda;
mu = p.PDE.mu;
x = p.level(end).x;

nrNodes = p.level(end).nrNodes;
nrEdges = p.level(end).nrEdges;
nrElems = p.level(end).nrElems;
gradC4e = p.level(end).enum.grad4e;
gradNC4e = p.level(end).enum.gradNC4e;
ed4e = p.level(end).enum.ed4e;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = x(1:nrNodes+nrEdges);
uC = u(1:nrNodes);
uNC = u(nrNodes+1:end);

gradU = zeros(2,2,nrElems);
epsilon = zeros(2,2,nrElems);
sigma = zeros(2,2,nrElems);

genericUnit = eye(2);

u4e = zeros(3,2,nrElems);

for curElem = 1:nrElems
    curNodes = n4e(curElem,:);
    curEdges = ed4e(curElem,:);
    gradC = gradC4e(:,:,curElem);
    gradNC = gradNC4e(:,:,curElem);

    u1 = uC(curNodes)';
    u2 = uNC(curEdges)';
     
	u4e(:,:,curElem) = [u1',u2'];
	
    uGrad1 = u1*gradC;
    uGrad2 = u2*gradNC;

    % calculation of piecewise constant gradient of u
    curGrad = [uGrad1;uGrad2];
    gradU(:,:,curElem) = curGrad;
	
    % calculation of piecewise constant strain (\eps(u))
    curStrain = (curGrad+curGrad')/2;
    epsilon(:,:,curElem) = curStrain;
    
    % calculation of piecewise constant stress (\sigma(u))
    traceUGrad = curGrad(1,1) + curGrad(2,2);
    curStress = lambda*traceUGrad*genericUnit + 2*mu*curStrain;
    sigma(:,:,curElem) = curStress;
end

%%%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%%%%%%%%%
p.level(end).u = u;
p.level(end).u4e = u4e;
p.level(end).gradU = gradU;
p.level(end).epsilon = epsilon;
p.level(end).sigma = sigma;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
