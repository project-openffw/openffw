function p = RT0P0postProc(p)
%postproc.m computes the post-processing datas for 
%the mixed RT0-P0-FE method.
%
%authors: David Guenther, Jan Reininghaus

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load geometry
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;

x = p.level(end).x;

% load enumerated data
ed4e = p.level(end).enum.ed4e;
e4ed = p.level(end).enum.e4ed;
length4ed = p.level(end).enum.length4ed;
area4e = p.level(end).enum.area4e;
nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = x(nrEdges+1:end)';
gradP = x(1:nrEdges);


grad = zeros(3,2,nrElems);
grad4e = zeros(nrElems,3);

for curElem = 1:nrElems
	curEdges = ed4e(curElem,:);
	curP = x(curEdges);
	curCoords = c4n(n4e(curElem,:),:)';
	area = area4e(curElem);
	
    grad4e(curElem,:) = curP;
    
	signum = ones(3,1);
	I = find(e4ed(curEdges,2) == curElem);
	signum(I) = -1;
		
	psi1P1 = signum(1)*length4ed(curEdges(1))/(2*area)*(curCoords(:,1)-curCoords(:,3));	
	psi3P1 = signum(3)*length4ed(curEdges(3))/(2*area)*(curCoords(:,1)-curCoords(:,2));
	
	psi1P2 = signum(1)*length4ed(curEdges(1))/(2*area)*(curCoords(:,2)-curCoords(:,3));	
	psi2P2 = signum(2)*length4ed(curEdges(2))/(2*area)*(curCoords(:,2)-curCoords(:,1));
	
	psi2P3 = signum(2)*length4ed(curEdges(2))/(2*area)*(curCoords(:,3)-curCoords(:,1));	
	psi3P3 = signum(3)*length4ed(curEdges(3))/(2*area)*(curCoords(:,3)-curCoords(:,2));
	

	grad(1,:,curElem) = (curP(1)*psi1P1 + curP(3)*psi3P1);
	grad(2,:,curElem) = (curP(1)*psi1P2 + curP(2)*psi2P2);
	grad(3,:,curElem) = (curP(2)*psi2P3 + curP(3)*psi3P3);
	
end

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).u = u;
p.level(end).u4e = u;
p.level(end).grad = grad;
p.level(end).grad4e = grad4e;
p.level(end).gradP = gradP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
