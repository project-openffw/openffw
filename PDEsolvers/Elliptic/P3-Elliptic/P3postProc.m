function p = P3postProc(p)

% author: Joscha Gedicke

%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = p.level(end).x;

% load enumerated data
nrElems = p.level(end).nrElems;
dofU4e = p.level(end).enum.dofU4e;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = x(:);

u4e = zeros(nrElems,size(dofU4e,2));

for curElem = 1:nrElems
	curNodes = dofU4e(curElem,:);
	curU = u(curNodes);
    u4e(curElem,:) = curU;
end

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).u = u;
p.level(end).u4e = u4e;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%