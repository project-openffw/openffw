function p = ODRTestimateDWRNEW(p)

% Copyright 2008 Joscha Gedicke, Lena Noack
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
% length4ed = p.level(end).enum.length4ed;
n4e = p.level(end).geom.n4e;
% n4ed = p.level(end).enum.n4ed;
% ed4e = p.level(end).enum.ed4e;
curLvl = length(p.level);
degree = 19;
nrNodes = p.level(end).nrNodes;


%% ESTIMATE
% h_T = max(length4ed(ed4e),[],2);

nu_T = integrate(n4e,curLvl,2*degree,@integrand,p);
primal= nu_T;

% weight
dnu_T = integrate(n4e,curLvl,2*degree,@DWRAverageError,p);
dual = dnu_T;


nu = primal.*dual;


%% OUTPUT
p.level(end).etaT = nu;
p.level(end).estimatedError = sum(sqrt(nu));

function val = integrand(x,y,parts,lvl,p)
gradUh = permute(p.statics.sigma_h(x,y,parts,lvl,p),[3 2 1]);
q = -[x./2 y./2];
val(1,1,:) = (gradUh(:,1)-q(:,1)).^2+(gradUh(:,2)-q(:,2)).^2;

% Averaging of dual Solution: Ap_h-p_h
function val = DWRAverageError(x,y,parts,lvl,p)
DWRAp_h = p.statics.Ap_h;
DWRp_h = p.statics.grad_h;
%DWRAp_h = p.statics.DWRAp_h;
%DWRp_h = p.statics.DWRgrad_h;

%% Residuum
DWRAph = DWRAp_h(x,y,parts,lvl,p);
DWRph = DWRp_h(x,y,parts,lvl,p);

residuum = DWRAph - DWRph;
%residuum = permute(residuum,[3 1 2]);

%% RETURN
val(1,:,:) = (sum(residuum.^2,2))';