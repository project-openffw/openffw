function p = P1CRestimate(p)
%estimate.m estimate the energy error for for the Kouhia-Stenberg FE 
%in linear elasticity through computing the residual and 
%the jump of the stress. 
%
%authors: David Guenther, Jan Reininghaus
%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Db = p.level(end).geom.Db;
Nb = p.level(end).geom.Nb;

f = p.problem.f;
g = p.problem.g;

uGrad = p.level(end).gradU;
sigma = p.level(end).sigma;

area4e = p.level(end).enum.area4e;
length4ed = p.level(end).enum.length4ed;
midPoint4e = p.level(end).enum.midPoint4e;
midPoint4ed = p.level(end).enum.midPoint4ed;
ed4e = p.level(end).enum.ed4e;
DbEd = p.level(end).enum.DbEd;
NbEd = p.level(end).enum.NbEd;
normals4e = p.level(end).enum.normals4e;
normals4NbEd = p.level(end).enum.normals4NbEd;
nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

etaT = zeros(nrElems,1);
etaEdNormal = zeros(nrEdges,2);
etaEdTangent = zeros(nrEdges,2);

f4e = sum( f(midPoint4e(:,1),midPoint4e(:,2),p).^2,2 );

for curElem = 1:nrElems
    curEdges = ed4e(curElem,:);
    curLengthEdge = length4ed(curEdges);
    longestEdge = max(curLengthEdge);
    area = area4e(curElem);
    normal = normals4e(:,:,curElem);
    tangent = normal(:,[2 1]).*[1 -1;1 -1;1 -1];
    curGrad = uGrad(:,:,curElem);
    curSigma = sigma(:,:,curElem);

    %%% ETA_T %%%%%% 1-point gauss quadrature %%%%%%%%%%%
    % eta_T^2 = h_T^2 * norm(f)^2_L2(T) = h_T^2 * f(x_T)^2 * area
    etaT(curElem) = longestEdge^2 * f4e(curElem) * area;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% ETA_Ed (normal)%%%%%% 1-point gauss quadrature %%%%%%%%%%%
    % eta_Ed^2 = h_E * norm(J(sig))^2_L2(E) = h_E^2 * (sig*n_1 + sig*n_2)^2
    etaEdNormal(curEdges,:) = etaEdNormal(curEdges,:) + (curSigma*normal')'.*([curLengthEdge,curLengthEdge]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% ETA_Ed (tangent)%%%%%% 1-point gauss quadrature %%%%%%%%%%%
    % eta_Ed^2 = h_E * norm(du/ds)^2_L2(E) = h_E^2 * (U_grad*s_1 + U_grad*s_2)^2
    etaEdTangent(curEdges,:) = etaEdTangent(curEdges,:) + (curGrad*tangent')'.*([curLengthEdge,curLengthEdge]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Neumann condition
if ~isempty(Nb)
    g4Nb = g(midPoint4ed(NbEd,1),midPoint4ed(NbEd,2),normals4NbEd,p);
    etaEdNormal(NbEd,:) = etaEdNormal(NbEd,:) - g4Nb.*[length4ed(NbEd),length4ed(NbEd)]/2;
end

% Dirichlet condition
if ~isempty(Db)
    etaEdNormal(DbEd,:) = 0;
    etaEdTangent(DbEd,:) = 0;
end

etaJNormal = etaEdNormal(:,1).^2 + etaEdNormal(:,2).^2;
etaJTangent = etaEdTangent(:,1).^2 + etaEdTangent(:,2).^2;

% eta = sqrt( sum_T(eta_T^2) + sum_Ed(eta_Ed^2) )
estimatedError = sqrt( ( sum(etaT) + sum(etaJNormal) + sum(etaJTangent) ));
                       
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).etaT  = etaT;
p.level(end).etaEd  = etaJNormal + etaJTangent;
p.level(end).estimatedError = estimatedError;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
