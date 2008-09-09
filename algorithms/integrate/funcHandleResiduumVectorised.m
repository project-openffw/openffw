function val = funcHandleResiduumVectorised(x,y,x_ref,y_ref,parts,curLvl,p)
% function handle for calculation of the residuum

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


%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kappa = p.problem.kappa;
lambda = p.problem.lambda;
mu = p.problem.mu;
f  = p.problem.f;
u_h    = p.statics.u_hVectorised;
gradU_h    = p.statics.gradU_hVectorised;
D2U_h    = p.statics.d2U_hVectorised;

%% assume piecewise constant coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
midPoint4e = p.level(curLvl).enum.midPoint4e(parts,:);
kappa4e = permute(kappa(midPoint4e(:,1),midPoint4e(:,2),p),[3 1 2]);
lambda4e = lambda(midPoint4e(:,1),midPoint4e(:,2),p);
mu4e = mu(midPoint4e(:,1),midPoint4e(:,2),p);

%% Residuum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curf = f(x,y,p);
curU_h = u_h(x,y,x_ref,y_ref,parts,curLvl,p);
curGradU_h = permute(gradU_h(x,y,x_ref,y_ref,parts,curLvl,p),[3 2 1]);
curD2U_h = permute(D2U_h(x,y,x_ref,y_ref,parts,curLvl,p),[3 1 2]);


%-div(kappa gradU_h) + lambda*gradU_h + mu u_h - eigenvalue_h*u_h
residuum = -kappa4e(:,1,1).*curD2U_h(:,1,1) - kappa4e(:,1,2).*curD2U_h(:,2,1) ...
           -kappa4e(:,2,1).*curD2U_h(:,1,2) - kappa4e(:,2,2).*curD2U_h(:,2,2) ...
           + lambda4e(:,1).*curGradU_h(:,1) + lambda4e(:,2).*curGradU_h(:,2) ...
           + mu4e.*curU_h(:) ...
           - curf(:);

%% Return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
val(1,1,:) = (residuum.^2)';
