function p = AWenumerate(p)
%enumerate.m creates all necessarily data 
%for the Arnold-Winther mixed FE in linear elasticity. 
%
%authors: David Guenther, Jan Reininghaus

p = genericEnumerate(p);

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;
nrNodes = p.level(end).nrNodes;
nrEdges = p.level(end).nrEdges;
nrElems = p.level(end).nrElems;
area4e = p.level(end).enum.area4e;
ed4e = p.level(end).enum.ed4e;
lambda = p.PDE.lambda;
mu = p.PDE.mu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gradients of P1-Hat functions
grad4e = getGrad4e(c4n,n4e,area4e);
nrDoF = 3*nrNodes+4*nrEdges+3*nrElems+6*nrElems;

massMatrix = 1/2520 * ...
 [  420    210    210     84     42     84     14      0    -14     14      
    210    420    210     84     84     42    -14     14      0     14      
    210    210    420     42     84     84      0    -14     14     14      
     84     84     42     28     14     14      0      2     -2      4      
     42     84     84     14     28     14     -2      0      2      4      
     84     42     84     14     14     28      2     -2      0      4      
     14    -14      0      0     -2      2      3     -1     -1      0      
      0     14    -14      2      0     -2     -1      3     -1      0      
    -14      0     14     -2      2      0     -1     -1      3      0      
     14     14     14      4      4      4      0      0      0      1  ];

C = mu*[2,0,0;0,2,0;0,0,1] + lambda*[1,1,0;1,1,0;0,0,0];

lame1 = (lambda + 2*mu)/(4*mu*(lambda+mu));
lame2 = -lambda /(4*mu*(lambda+mu));
lame3 = 1/(2*mu);

Cinv = [lame1 lame2 0;
        lame2 lame1 0
         0      0   lame3];
     
dofU4e = zeros(nrElems,6);
dofSigma4e = zeros(nrElems,24);  

for curElem = 1:nrElems
    curNodes = n4e(curElem,:);
    curEdges = ed4e(curElem,:);
    
    dof4Nodes = repmat(3*curNodes,3,1) - repmat([2;1;0],1,3);
    dof4Edges = 3*nrNodes + (repmat(4*curEdges,4,1) - repmat([3;2;1;0],1,3));
    dof4Elem  = 3*nrNodes + 4*nrEdges + 3*(curElem-1) + (1:3);
    dof4Displ = 3*nrNodes + 4*nrEdges + 3*nrElems + 6*(curElem-1) + (1:6);

    dof4Nodes = dof4Nodes(:);
    dof4Edges = dof4Edges(:);
    
    dofU4e(curElem,:) = dof4Displ';
    dofSigma4e(curElem,:) = [dof4Nodes; dof4Edges; dof4Elem']'; 
end

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).enum.dofU4e = dofU4e;
p.level(end).enum.dofSigma4e = dofSigma4e;
p.level(end).enum.grad4e = grad4e;
p.level(end).enum.massMatrix = massMatrix;
p.level(end).enum.C = C;
p.level(end).enum.Cinv = Cinv;
p.level(end).nrDoF = nrDoF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = AWgetBasisCoefficents(p);
