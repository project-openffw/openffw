function p = TWRTprolong(p,lvl)
% author: David Guenther 

%% Input
parents4e = p.level(lvl).enum.parents4e;
e4ed = p.level(lvl).enum.e4ed;
ed4e = p.level(lvl).enum.ed4e;
normal4ed = p.level(lvl).enum.normals4ed;
normals4e = p.level(lvl).enum.normals4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;
midPoint4e = p.level(lvl).enum.midPoint4e;
nrEdges = p.level(lvl).nrEdges;
nrElems = p.level(lvl).nrElems;

sigma_h = p.statics.sigma_h;
DWRsigma_h = p.statics.DWRsigma_h;
u_h = p.statics.u_h;
DWRu_h = p.statics.DWRu_h;

x0 = loadField('p.epsilon(end).','x0',p);
DWRx0 = loadField('p.epsilon(end).','DWRx0',p);

%% prolongation

if lvl == 2
    x0 = ones(nrEdges+nrElems,1);
    DWRx0 = ones(nrEdges+nrElems,1);
    p.level(lvl).x = x0;
    p.level(lvl).DWRx = DWRx0;
    return
end

if isempty(x0)
x0 = ones(nrEdges+nrElems,1);
DWRx0 = ones(nrEdges+nrElems,1);

    if lvl > 2
        % RT0-prolongation
        for curEdge = 1:nrEdges
            elems = e4ed(curEdge,:);
            if elems(2) == 0
                elems(2) = elems(1);
            end
            parents = parents4e(elems);

            midPoint = midPoint4ed(curEdge,:);
            normal = normal4ed(curEdge,:);

            %%%%
            edges = ed4e(elems,:);
            index1 = find(curEdge == edges(1,:));
            index2 = find(curEdge == edges(2,:));
            normal1 = normals4e(index1,:,elems(1));
            normal2 = normals4e(index2,:,elems(2));
            %%%%%

            sigma4T1 = sigma_h(midPoint(1),midPoint(2),parents(1),lvl-1,p)*normal1';
            sigma4T2 = sigma_h(midPoint(1),midPoint(2),parents(2),lvl-1,p)*normal1';
            DWRsigma4T1 = DWRsigma_h(midPoint(1),midPoint(2),parents(1),lvl-1,p)*normal1';
            DWRsigma4T2 = DWRsigma_h(midPoint(1),midPoint(2),parents(2),lvl-1,p)*normal1';

            DWRx0(curEdge) = 1/2*(DWRsigma4T1 + DWRsigma4T2);
        end


        %P0-prolongation
        for curElem = 1:nrElems
            parent = parents4e(curElem);
            midPoint = midPoint4e(curElem,:);
            curU = u_h(midPoint(1),midPoint(2),parent,lvl-1,p);
            DWRcurU = DWRu_h(midPoint(1),midPoint(2),parent,lvl-1,p);
            x0(nrEdges+curElem) = curU;
            DWRx0(nrEdges+curElem) = DWRcurU;
        end
    end

end
%% OUTPUT
p.level(lvl).x = x0;
p.level(lvl).DWRx = DWRx0;
