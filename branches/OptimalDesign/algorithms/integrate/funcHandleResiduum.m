function val = funcHandleResiduum(x,y,curElem,curLvl,p)

% author: Joscha Gedicke

%% INPUT
f = p.problem.f;
kappa = p.problem.kappa;
lambda = p.problem.lambda;
mu = p.problem.mu;
u_h    = p.statics.u_h;
gradU_h    = p.statics.gradU_h;
D2U_h    = p.statics.D2U_h;

%% assume piecewise constant coefficients
midPoint4e = p.level(curLvl).enum.midPoint4e;
curMidPoint = midPoint4e(curElem,:);
curKappa = kappa(curMidPoint(1),curMidPoint(2),p);
curLambda = lambda(curMidPoint(1),curMidPoint(2),p);
curMu = mu(curMidPoint(1),curMidPoint(2),p);

%% Residuum
curU_h = u_h(x,y,curElem,curLvl,p);
curGradU_h = gradU_h(x,y,curElem,curLvl,p);
curD2U_h = D2U_h(x,y,curElem,curLvl,p);
curf     = f(x,y,p);

%-div(kappa gradU_h) + lambda*gradU_h + mu u_h - f
residuum = -curKappa(1,1)*curD2U_h(:,1,1) - curKappa(1,2)*curD2U_h(:,1,2) ...
           -curKappa(2,1)*curD2U_h(:,2,1) - curKappa(2,2)*curD2U_h(:,2,2) ...
           + curLambda(1)*curGradU_h(:,1) + curLambda(2)*curGradU_h(:,2) ...
           + curMu*curU_h(:) ...
           - curf(:);

%% RETURN
val(1,:,:) = (residuum.^2)';
