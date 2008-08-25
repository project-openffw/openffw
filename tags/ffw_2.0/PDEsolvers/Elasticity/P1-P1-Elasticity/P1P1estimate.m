function p = P1P1estimate(p)
% estimate the energy error for a P1-FE method 
% in linear elasticity through computing the residual and 
% the jump of the stress. 

% Copyright 2007 Jan Reininghaus, David Guenther
%
% This file is part of FFW.
%
% FFW is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% FFW is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Db = p.level(end).geom.Db;
Nb = p.level(end).geom.Nb;

f = p.problem.f;
g = p.problem.g;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta_T = zeros(nrElems,1);
etaEdNormal = zeros(nrEdges,2);

f4e = sum( f(midPoint4e,p).^2,2 );

for curElem = 1:nrElems
    curEdges = ed4e(curElem,:);
    curLengthEdge = length4ed(curEdges);
    longestEdge = max(curLengthEdge);
    area = area4e(curElem);
    curSigma = sigma(:,:,curElem);
    
    normal = normals4e(:,:,curElem);
    
    %%% ETA_T %%%%%% 1-point gauss quadrature %%%%%%%%%%%
    % eta_T^2 = h_T^2 * norm(f)^2_L2(T) = h_T^2 * f(x_T)^2 * area
    eta_T(curElem) = longestEdge^2 * f4e(curElem) * area;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% ETA_Ed %%%%%% 1-point gauss quadrature %%%%%%%%%%%
    % eta_Ed^2 = h_E * norm(J(sig))^2_L2(E) = h_E^2 * (sig*n_1 + sig*n_2)^2
    etaEdNormal(curEdges,:) = etaEdNormal(curEdges,:) + (curSigma*normal')'.*([curLengthEdge,curLengthEdge]);   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Neumann condition
if ~isempty(Nb)
    g4Nb = g(midPoint4ed(NbEd,:),normals4NbEd,p);
    etaEdNormal(NbEd,:) = etaEdNormal(NbEd,:) - g4Nb.*[length4ed(NbEd),length4ed(NbEd)]/2;
end

% Dirichlet condition
if ~isempty(Db)
    etaEdNormal(DbEd,:) = 0;
end

% J(sig)^2 = J(sig)*J(sig) (scalar product)
eta_Jn = etaEdNormal(:,1).^2 + etaEdNormal(:,2).^2;
    
% eta = sqrt( sum_T(eta_T^2) + sum_Ed(eta_Ed^2) )
estimatedError = sqrt( sum(eta_T) + sum(eta_Jn) );    

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).etaT  = eta_T;
p.level(end).etaEd  = eta_Jn;
p.level(end).estimatedError = estimatedError;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
