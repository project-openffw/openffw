function p = P1postProc(p)
%postproc.m computes the post-processing datas for 
%a conforming P1-FE method.
%
%authors: David Guenther, Jan Reininghaus
%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load geometry
n4e = p.level(end).geom.n4e;

x = p.level(end).x;

% load enumerated data
P1grad4e = p.level(end).enum.grad4e;
nrNodes = p.level(end).nrNodes;
nrElems = p.level(end).nrElems;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = x(1:nrNodes);

grad4e = zeros(nrElems,2);
u4e = zeros(nrElems,3);

for curElem = 1:nrElems
	curGrads = P1grad4e(:,:,curElem);
	curNodes = n4e(curElem,:);
	curU = u(curNodes);
    u4e(curElem,:) = curU;
	grad4e(curElem,:) = sum([curU,curU] .* curGrads,1);
end

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).u = u;
p.level(end).u4e = u4e;
p.level(end).grad4e = grad4e;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%