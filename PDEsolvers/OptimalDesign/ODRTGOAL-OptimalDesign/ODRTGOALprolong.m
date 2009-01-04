function p = ODRTGOALprolong(p,lvl)
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
lambda1 = p.statics.lambda1;
lambda2 = p.statics.lambda2;
u_h = p.statics.u_h;

%% prolongation
x0 = zeros(2*nrEdges+2*nrElems,1);

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

        sigma4T1 = sigma_h(midPoint,parents(1),lvl-1,p)*normal1';
        sigma4T2 = sigma_h(midPoint,parents(2),lvl-1,p)*normal1';
        
        lambda14T1 = lambda1(midPoint,parents(1),lvl-1,p)*normal1';
        lambda14T2 = lambda1(midPoint,parents(2),lvl-1,p)*normal1';
        
        x0(curEdge) = 1/2*(sigma4T1 + sigma4T2);
        x0(nrEdges+nrElems+curEdge) = 1/2*(lambda14T1 + lambda14T2);
    end


    %P0-prolongation
    for curElem = 1:nrElems
        parent = parents4e(curElem);
        midPoint = midPoint4e(curElem,:);
        curU = u_h(midPoint,parent,lvl-1,p);
        curLambda = lambda2(midPoint,parent,lvl-1,p);
        x0(nrEdges+curElem) = curU;
        x0(nrEdges+nrElems+nrEdges) = curLambda;
    end
end

%% OUTPUT
p.level(lvl).x = x0;
