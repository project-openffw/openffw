function val = funcHandleResiduum(x,y,curElem,curLvl,p)
% function handle to calculate the residuum

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


%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = p.problem.f;
kappa = p.problem.kappa;
lambda = p.problem.lambda;
mu = p.problem.mu;
u_h    = p.statics.u_h;
gradU_h    = p.statics.gradU_h;
D2U_h    = p.statics.d2U_h;

%% assume piecewise constant coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
midPoint4e = p.level(curLvl).enum.midPoint4e;
curMidPoint = midPoint4e(curElem,:);
curKappa = kappa(curMidPoint(1),curMidPoint(2),p);
curLambda = lambda(curMidPoint(1),curMidPoint(2),p);
curMu = mu(curMidPoint(1),curMidPoint(2),p);

%% Residuum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curU_h = u_h(x,y,curElem,curLvl,p);
curGradU_h = permute(gradU_h(x,y,curElem,curLvl,p),[3 2 1]);
curD2U_h = permute(D2U_h(x,y,curElem,curLvl,p),[3 1 2]);
curf     = f(x,y,p);

%-div(kappa gradU_h) + lambda*gradU_h + mu u_h - f
residuum = -curKappa(1,1)*curD2U_h(:,1,1) - curKappa(1,2)*curD2U_h(:,2,1) ...
           -curKappa(2,1)*curD2U_h(:,1,2) - curKappa(2,2)*curD2U_h(:,2,2) ...
           + curLambda(1)*curGradU_h(:,1) + curLambda(2)*curGradU_h(:,2) ...
           + curMu*curU_h(:) ...
           - curf(:);

%% Return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
val(1,:,:) = (residuum.^2)';
