function val = DWRDualAverage(pts,pts_ref,parts,curLvl,p)

% Copyright 2007 Joscha Gedicke, Lena Noack
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
fluxRecovery   = p.statics.fluxRecoveryUh;
sigma_h = p.statics.sigma_h;

%% Residuum
fluxRecovery  = fluxRecovery(pts,pts_ref,parts,curLvl,p);
sigma_h  = sigma_h(pts,pts_ref,parts,curLvl,p);
sigma_h  = permute(sigma_h,[3 2 1]);

%tau_h - A tau_h
sigma_h = reshape(sigma_h, [2 1 size(sigma_h,3)]);
residuum = fluxRecovery-sigma_h;
residuum = permute(residuum,[3 1 2]);

%% RETURN
val(1,:,:) = (sum(residuum.^2,2))';

