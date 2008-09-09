function p = P1P0createLinSys(p)
% just pure natural Dirichlet boundary condition
% P0-P1 method is unstable

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load geometry
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;

Nb = p.level(end).geom.Nb;

% load problem definition
f = p.problem.f;
g = p.problem.g;
u_D = p.problem.u_D;

% load enumerated data

midPoint4ed = p.level(end).enum.midPoint4ed;
midPoint4e = p.level(end).enum.midPoint4e;
length4ed = p.level(end).enum.length4ed;
e4ed = p.level(end).enum.e4ed;
DbEd = p.level(end).enum.DbEd;
NbEd = p.level(end).enum.NbEd;
ed4e = p.level(end).enum.ed4e;
area4e = p.level(end).enum.area4e;
normals4NbEd = p.level(end).enum.normals4NbEd;

nrNodes = p.level(end).nrNodes;
nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assembling global energy matrix				   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = zeros(6,6,nrElems);
C = zeros(1,6,nrElems);

genericMass = 1/12* [   2 1 1 0 0 0;
                        1 2 1 0 0 0;
                        1 1 2 0 0 0;
                        0 0 0 2 1 1;
                        0 0 0 1 2 1;
                        0 0 0 1 1 2];

% genericC = 2*[1 0 -1 0 1 -1];
genericC = [1 0 -1 0 1 -1];

for curElem = 1:nrElems
% 	curNodes = n4e(curElem,:);
% 	curCoords = c4n(curNodes([3 1 2]),:)';
	area = area4e(curElem);
	B(:,:,curElem) = area*genericMass;
    C(:,:,curElem) = area*genericC;
end

% local DoF -> global DoF marix B
globalDoF = [n4e,nrNodes + n4e];
globalDoF_T = globalDoF';
I = repmat(globalDoF_T(:),1,size(globalDoF,2))';
J = repmat(globalDoF,1,size(globalDoF,2))';

B = sparse(I(:),J(:),B(:));

I = repmat((1:nrElems),6,1);
J = [n4e,nrNodes + n4e]';
C = sparse(I(:),J(:),C(:),nrElems,2*nrNodes);

A = [B		C';
	 C		sparse(nrElems,nrElems)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assembling Righthandside					   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = zeros(2*nrNodes + nrElems,1);
b(2*nrNodes+1:end) = -area4e .* f(midPoint4e(:,1),...
								midPoint4e(:,2),p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Include boundary conditions						  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = zeros(2*nrNodes + nrElems,1);

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).A = A;
p.level(end).b = b;
p.level(end).x = x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
