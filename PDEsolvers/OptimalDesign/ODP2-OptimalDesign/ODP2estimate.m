function p = ODP2estimate(p)
% author: David Guenther 
% Copyright 2007 David Guenther
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
%
%
%

%% Input
length4ed = p.level(end).enum.length4ed;
n4ed = p.level(end).enum.n4ed;
degree = p.params.rhsIntegtrateExactDegree;

curLvl = size(p.level,2);

%% compute jump error and oscillations
eta4ed = length4ed.*integrate(n4ed,curLvl,degree,@jumpError,p);

%% compose error terms
estimatedError = sqrt( sum(eta4ed) );

%% Output
p.level(end).etaEd  = sqrt(eta4ed);
p.level(end).estimatedError = estimatedError;

%% supply jump error
function val = jumpError(points,curEdge,lvl,p)

sigma_h = p.statics.sigma_h;
e4ed = p.level(lvl).enum.e4ed;
normals4ed = p.level(lvl).enum.normals4ed;

elems = e4ed(curEdge,:);
normal = normals4ed(curEdge,:);

evalSigma1 = sigma_h(points,elems(1),lvl,p)*normal';

if elems(2) ~= 0
    % E is an interior edge
    evalSigma2 = sigma_h(points,elems(2),lvl,p)*normal';
else
    % E is a Dirichlet Edge
    evalSigma2 = evalSigma1;
end

val = zeros(1,1,length(points(:,1)));
val(1,1,:) = sum((evalSigma1 - evalSigma2).^2,2);
