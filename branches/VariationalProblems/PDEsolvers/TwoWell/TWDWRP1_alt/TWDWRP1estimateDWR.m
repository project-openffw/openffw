function p = TWDWRP1estimateDWR_Jump_SigmaP(p)
% Copyright 2007 Joscha Gedicke, Lena NOack

%% INPUT
length4ed = p.level(end).enum.length4ed;
n4e = p.level(end).geom.n4e;
n4ed = p.level(end).enum.n4ed;
ed4e = p.level(end).enum.ed4e;
curLvl = length(p.level);
degree = 1;


%% ESTIMATE
nu_T = 1;%integrateVectorised(n4e,curLvl,2*degree,@funcHandleEigenvalueResiduumVectorised,p);
nu_E = 1,%integrateVectorised(n4ed,curLvl,2*degree,@funcHandleNormalJumpVectorised,p);
% primalResiduum = sqrt( nu_T + 1/2./h_T.*sum(nu_E(ed4e),2) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F1 = -(1/sqrt(13))*[3;2];
F2 = (1/sqrt(13))*[3;2];
options = optimset('Display','off','Jacobian','on','NonlEqnAlgorithm','gn','MaxIter',10,'TolFun',1e-10,'TolX',1e-10);
p2 = initFFW('P1',p.params.problem.name,'bulk',100,'TwoWell','TWdirectFsolve','redGreenBlue','Dual');
p.params.options = options;
p.params.RHS = 'RHS2';
p.params.CONV = 'c'; %convex energy density
%p.params.CONV = 'nc'; %non-convex energy density
p.problem.F1 = F1;
p.problem.F2 = F2;
p.params.rhsIntegtrateExactDegree = 19;
p.params.nonLinearExactIntegrateDegree = 5;
p = p.statics.run(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p2 = initFFW('P1',p.params.problem.name,'bulk',100,'TwoWell','TWdirectFsolve','redGreenBlue','Dual');
%p2 = initFFW('P1',p.params.problem.name,'P1Dual',inf,'TwoWell',p.params.solver,'P1Dual','Dual');
p2.level(1).geom = p.level(curLvl).geom;
p2.params.modules.mark.refineFirstLevels = 0;
p2.params.maxLevel = 1;
p2 = computeSolution(p2);
p.level(curLvl).P2u4e = p2.level(end).u4e;
% p.level(curLvl).u4e = p2.level(end).u4e(:,1:3);
mu_T = 1;%integrateVectorised(n4e,curLvl,4*degree,@funcHandleEigenvalueWeightElemsVectorised,p);
mu_E = 1;%sintegrateVectorised(n4ed,curLvl,4*degree,@funcHandleEigenvalueWeightEdgesVectorised,p);
% primalWeight = sqrt( mu_T + 1/2*h_T.*sum(mu_E(ed4e),2) );

% dual residuum
q = p;
q.level(end).A = q.level(end).A';
q = q.statics.solve(q);
q = q.statics.postProc(q); 
dnu_T = 1;%integrateVectorised(n4e,curLvl,2*degree,@funcHandleEigenvalueResiduumVectorised,q);
dnu_E = 1;%integrateVectorised(n4ed,curLvl,2*degree,@funcHandleNormalJumpVectorised,q);
% dualResiduum = sqrt( dnu_T + 1/2./h_T.*sum(dnu_E(ed4e),2) );

q2 = p2;
q2.level(end).A = q2.level(end).A';
q2 = q2.statics.solve(q2);
q2 = q2.statics.postProc(q2); 
q.level(curLvl).P2u4e = q2.level(end).u4e;
% q.level(curLvl).u4e = q2.level(end).u4e(:,1:3);
dmu_T = 1;%integrateVectorised(n4e,curLvl,4*degree,@funcHandleEigenvalueWeightElemsVectorised,q);
dmu_E = 1;%integrateVectorised(n4ed,curLvl,4*degree,@funcHandleEigenvalueWeightEdgesVectorised,q);
% dualWeight = sqrt( dmu_T + 1/2*h_T.*sum(dmu_E(ed4e),2) );

% nu = primalResiduum.*dualWeight+dualResiduum.*primalWeight;
nu = sqrt(nu_T).*sqrt(dmu_T)+1/2.*sqrt(sum(nu_E(ed4e),2)).*sqrt(sum(dmu_E(ed4e),2))...
    +sqrt(dnu_T).*sqrt(mu_T)+1/2.*sqrt(sum(dnu_E(ed4e),2)).*sqrt(sum(mu_E(ed4e),2));

%% OUTPUT
p.level(end).etaT = nu;
p.level(end).estimatedError = sum(nu);

