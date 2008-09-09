function p = AWpostProc(p)
%postproc.m computes the post-processing datas 
%for the Arnold-Winther mixed FE in linear elasticity. 
%
%authors: David Guenther, Jan Reininghaus
%%%%%%%%%%%%%%  INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = p.PDE.lambda;
mu = p.PDE.mu;

x = p.level(end).x;
A = p.level(end).A;

nrNodes = p.level(end).nrNodes;
nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;
nrDoF = p.level(end).nrDoF;
area4e = p.level(end).enum.area4e;
dofSigma4e = p.level(end).enum.dofSigma4e;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


index1 = 3*nrNodes + 4*nrEdges + 3*nrElems + 1;
index2 = 3*nrNodes + 4*nrEdges + 3*nrElems + 6*nrElems;

dummy = repmat(sqrt(area4e)',6,1);
nrSigmaDoF = 3*nrNodes + 4*nrEdges + 3*nrElems;
displacementScaling = [ones(nrSigmaDoF,1);...
                        dummy(:);...
                        ones(size(A,1) - nrSigmaDoF - 6*nrElems,1)];
S = spdiags(displacementScaling,0,size(A,1),size(A,2));
x = x*S;

% get discrete sigma
sigma = zeros(nrElems,24);
for curElem = 1:nrElems	
    curSigmaDoFs = dofSigma4e(curElem,:);
    sigma(curElem,:) = x(curSigmaDoFs)';
end

%get the energy of the discrete problem
dof4p = 3*nrNodes+4*nrEdges+3*nrElems;
J = 1:dof4p;
energy = sqrt(  x(J)*A(J,J)*x(J)'   );

%compute deviatoric part of u
dev = zeros(nrNodes,1);
for j = 1:nrNodes
    dev(j) = ( (mu^2/(6*(mu+lambda)^2)+0.5)*(x(3*j-2)+x(3*j-1))^2+...
                2*(x(3*j)^2-x(3*j-2)*x(3*j-1)) )/4/mu;
end

%get the displacement u of each node of a triangle 
u = x(dof4p+1:dof4p+6*nrElems);
u = reshape(u,3,2,nrElems); 

%%%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).x = x;
p.level(end).u4e = u;
p.level(end).sigma = sigma;
p.level(end).dev4n = dev;
p.level(end).energy = energy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
