function p = ODRTestimate_Proj(p)
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

%% INPUT
n4e = p.level(end).geom.n4e;
lvl = size(p.level,2);
degree = loadField('p.params','nonLinearExactIntegrateDegree',p,19);

%% compute (\grad v_h,\grad w_h)-projection of DW*(p_h)
% get rhs for projection
f4e = integrate(n4e,lvl,degree,@integrandProj,p);
f4e = squeeze(f4e);
% f4ed = integrate(n4ed,lvl,degree,@integrandProjEdges,p);
% f4ed = squeeze(f4ed);
    
% set up projection method ('P1','P2')
projMethod4estimate = loadField('p.params','projMethod4estimate',p,'P1');
% initialize projection method
q = initFFW(projMethod4estimate,p.params.problem.name,'uniform',p.level(end).nrDoF,'optimalDesign','direct','redGreenBlue','STUBestimate');
q.problem.kappa_dummy = @(x,y,p) reshape([1+x-x, x-x, y-y 1+y-y]',[2,2,length(x)]);
q.problem.kappa = @kappa;
q.problem.lambda_dummy = @(x,y,p) [x-x, y-y];
q.problem.lambda = @lambda;
q.problem.mu_dummy = @(x,y,p) x-x;
q.problem.mu = @mu;
q.problem.g = [];
q.level(1).geom = p.level(end).geom;
q.level(1).f4e = f4e;
% q.level(1).f4ed = f4ed;
q.params.maxLevel = 1;
q.params.maxNrDoF = p.level(lvl).nrDoF+1;
q.params.modules.mark.refineFirstLevels = false;
%compute projection
q = elliptic_init(q);
q = computeSolution(q);
p.level(end).projStruct = q;

%% estimate the error by computing \eta=||DW*(p_h)-\grad v_h||_L^2
eta4T = integrate(n4e,lvl,degree,@integrandError,p);

%% compose the error terms
estimatedError = sqrt( sum(eta4T) );

%% OUTPUT
p.level(end).etaT = sqrt(eta4T);
p.level(end).estimatedError = estimatedError;

%% supply the integrand: (DW*(p_h)-\grad v_h)^2
function val = integrandError(points,curElem,lvl,p)

q = p.level(lvl).projStruct;
gradProj = q.statics.gradU_h;
grad_h = p.statics.grad_h;
area4e = p.level(lvl).enum.area4e;

evalGrad = grad_h(points,curElem,lvl,p);
evalGradProj = squeeze(gradProj(points,curElem,1,q))';

val = sum( (evalGrad-evalGradProj).^2,2 );
val = reshape(val,[1 1 length(points(:,1))]);

%% supply the integrand: DW*(p_h)*\grad w_h)^2
function val = integrandProj(points,curElem,lvl,p)

projMethod4estimate = loadField('p.params','projMethod4estimate',p,'P1');
gradBasis = loadField('p.statics',strcat('gradBasis',projMethod4estimate),p,[]);
grad_h = p.statics.grad_h;

evalGrad = grad_h(points,curElem,lvl,p);
evalGradBasis = gradBasis(points,curElem,lvl,p);

evalGrad = reshape(evalGrad',[2 1 length(points(:,1))]);
val = matMul(evalGradBasis,evalGrad);

%% supply integrand edges: DW^*(p_h)*\nu * w_h
function val = integrandProjEdges(points,curEdge,lvl,p)

projMethod4estimate = loadField('p.params','projMethod4estimate',p,'P1');
basis = loadField('p.statics',strcat('basis',projMethod4estimate),p,[]);
grad_h = p.statics.grad_h;
normals4ed = p.level(lvl).enum.normals4ed;
e4ed = p.level(lvl).enum.e4ed;

elems = e4ed(curEdge,:);
if elems(2) == 0
    elems(2) = elems(1);
end
normal = normals4ed(curEdge,:);

evalGrad_1 = grad_h(points,elems(1),lvl,p);
evalBasis_1 = basis(points,elems(1),lvl,p);
evalGrad_2 = grad_h(points,elems(2),lvl,p);
evalBasis_2 = basis(points,elems(2),lvl,p);

nrBasis = size(evalBasis_1,2);

gradNormal_1 = evalGrad_1*normal';
gradNormal_2 = evalGrad_2*normal';

gradNormalBasis_1 = (gradNormal_1*ones(1,nrBasis)).*evalBasis_1;
gradNormalBasis_2 = (gradNormal_2*ones(1,nrBasis)).*evalBasis_2;

val = gradNormalBasis_1 - gradNormalBasis_2;

val = reshape(val',[nrBasis,1,length(points(:,1))]);

%% function wrapper for elliptic pde-coefficents
function val = kappa(points,p)

x = points(:,1);
y = points(:,2);

val = p.problem.kappa_dummy(x,y,p);

function val = lambda(points,p)

x = points(:,1);
y = points(:,2);

val = p.problem.lambda_dummy(x,y,p);

function val = mu(points,p)

x = points(:,1);
y = points(:,2);

val = p.problem.mu_dummy(x,y,p);
