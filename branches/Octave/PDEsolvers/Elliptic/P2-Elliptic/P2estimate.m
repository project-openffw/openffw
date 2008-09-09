function p = P2estimate(p)
% calculates lokal error indicators

% Copyright 2007 Joscha Gedicke
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


%% INPUT
length4ed = p.level(end).enum.length4ed;
n4e = p.level(end).geom.n4e;
n4ed = p.level(end).enum.n4ed;
ed4e = p.level(end).enum.ed4e;
curLvl = length(p.level);
degree = loadField('p.params','rhsIntegtrateExactDegree',p,2);


%% ESTIMATE
h_T = max(length4ed(ed4e),[],2);

nu_T = h_T.^2.*integrateVectorised(n4e,curLvl,max(4,degree),@funcHandleResiduumVectorised,p);
nu_E = length4ed.*integrateVectorised(n4ed,curLvl,2,@funcHandleNormalJumpVectorised,p);

nu = sqrt( nu_T + 1/2*sum(nu_E(ed4e),2) );

%% OUTPUT
p.level(end).etaT = nu;
p.level(end).estimatedError = norm(nu);
