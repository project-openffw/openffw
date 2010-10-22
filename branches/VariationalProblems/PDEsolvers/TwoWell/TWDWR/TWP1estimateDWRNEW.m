function p = TWP1estimateCompareSols(p)

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
length4ed = p.level(end).enum.length4ed;
n4e = p.level(end).geom.n4e;
% n4ed = p.level(end).enum.n4ed;
ed4e = p.level(end).enum.ed4e;
curLvl = length(p.level);
degree = 19;
nrNodes = p.level(end).nrNodes;
nrEdges = p.level(end).nrEdges;
nrElems = p.level(end).nrElems;


%% ESTIMATE
%h_T = max(length4ed(ed4e),[],2);

nu_T = integrate(n4e,curLvl,2*degree,@integrand,p);
primal= nu_T;

dnu_T = integrate(n4e,curLvl,2*degree,@Average,p);
dual = dnu_T;

nu = sqrt(primal.*dual);

estErr = integrate(n4e,curLvl,2*degree,@integrandCalc,p);
nu = abs(estErr);

%% OUTPUT
p.level(end).etaT = nu;
p.level(end).estimatedError = abs(sum(estErr));
%p.level(end).estimatedError = sum(nu);
%p.level(end).estimatedError = norm(nu,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = integrand(x,y,parts,lvl,p)
sigma_h = permute(p.statics.sigma_h(x,y,parts,lvl,p),[3 2 1]);

q = -[x./2 y./2];
val(1,1,:) = (sigma_h(:,1)-q(:,1)).^2+(sigma_h(:,2)-q(:,2)).^2;

% Averaging of dual Solution: Av_h-v_h
function val = Average(x,y,parts,lvl,p)
Ap_h = p.statics.DWRAGrad_h;
%Ap_h = p.statics.DWRAp_h;
p_h = p.statics.DWRgrad_h;

%% Residuum
Aph = Ap_h(x,y,parts,lvl,p);
ph = p_h(x,y,parts,lvl,p);

residuum = Aph - ph;
%residuum = permute(residuum,[3 1 2]);

%% RETURN
val(1,:,:) = (sum(residuum.^2,2))';




function val = integrandCalc(x,y,parts,lvl,p)
sigma_h = permute(p.statics.sigma_h(x,y,parts,lvl,p),[3 2 1]);

Ap_h = p.statics.DWRAGrad_h;
%Ap_h = p.statics.DWRAp_h;
p_h = p.statics.DWRgrad_h;

Aph = Ap_h(x,y,parts,lvl,p);
ph = p_h(x,y,parts,lvl,p);

q = -[x./2 y./2];

%res1 = sqrt((sigma_h(:,1)-q(:,1)).^2 + (sigma_h(:,2)-q(:,2)).^2);
%res2 = sqrt((Aph(:,1) - ph(:,1)).^2 + (Aph(:,2) - ph(:,2)).^2);
%val(1,1,:) = (sum(res1.*res2,2))';

sigma_h = reshape(sigma_h,[size(ph,1),2]);

res1 = q - sigma_h;
res2 = Aph - ph;
residual = res1(:,1).*res2(:,1) + res1(:,2).*res2(:,2);
val(1,:,:) = (sum(residual,2))';
%val(1,1,:) = ((sigma_h(:,1)-q(:,1)).*(Aph(:,1) - ph(:,1)) + (sigma_h(:,2)-q(:,2)).*(Aph(:,2) - ph(:,2)));
%val(1,1,:) = residual;