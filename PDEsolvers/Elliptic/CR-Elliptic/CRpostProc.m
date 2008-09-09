function p = CRpostProc(p)
%postproc.m computes the post-processing datas for 
%a nonconforming CR-FE method.
%
%authors: David Guenther, Jan Reininghaus

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load geometry
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;

x = p.level(end).x;
f = p.problem.f;

% load enumerated data
NCgrad4e = p.level(end).enum.gradNC4e;
ed4e = p.level(end).enum.ed4e;
midPoint4e = p.level(end).enum.midPoint4e;
nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = x(1:nrEdges);

grad4e = zeros(nrElems,2);
u4e = zeros(nrElems,3);

for curElem = 1:nrElems
	curGrads = NCgrad4e(:,:,curElem);
	curEdges = ed4e(curElem,:);
	curU = u(curEdges);
    u4e(curElem,:) = curU;
    grad4e(curElem,:) = curU' * curGrads;
end

gradRT0 = zeros(3,2,nrElems);
% assume f is piecewise constant
f4T = f(midPoint4e(:,1),midPoint4e(:,2),p);
for curElem = 1:nrElems
	curGrads = repmat(grad4e(curElem,:),3,1);
	curNodes = n4e(curElem,:);
	midPoint = repmat(midPoint4e(curElem,:),3,1);
	gradRT0(:,:,curElem) = curGrads - f4T(curElem)*(c4n(curNodes,:)-midPoint)/2;
end

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).u = u;
p.level(end).u4e = u4e;
p.level(end).grad4e = grad4e;
p.level(end).gradRT0 = gradRT0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%