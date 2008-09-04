function p = ODP2prolong(p,lvl)
% author: David Guenther 

%% INPUT
nrNodes = p.level(end).nrNodes;
nrEdges = p.level(end).nrEdges;

%% prolongation
n4e = p.level(end).geom.n4e;
lvl = size(p.level,2);
degree = loadField('p.params','nonLinearExactIntegrateDegree',p,19);

x0 = zeros(nrNodes+nrEdges,1);
if lvl > 2
    q = initFFW('P2',p.params.problem.name,'uniform',p.level(end).nrDoF,'elliptic','direct','redGreenBlue','STUBestimate');
    q.problem.kappa = @(x,y,p) reshape([1+x-x, x-x, y-y 1+y-y]',[2,2,length(x)]);
    q.problem.lambda = @(x,y,p) [x-x, y-y];
    q.problem.mu = @(x,y,p) x-x;
    q.problem.g = [];
    q.level(1).geom = p.level(end).geom;
    f4e = integrate(n4e,lvl,degree,@integrandProj,p);
    f4e = squeeze(f4e);
    q.level(1).f4e = f4e;
    q.params.maxLevel = 1;
    q.params.maxNrDoF = p.level(lvl).nrDoF+1;
    q.params.modules.mark.refineFirstLevel = false;
    %compute projection
    q = q.statics.run(q);
    x0 = q.level(1).x;
end

%% OUTPUT
p.level(lvl).x = x0;

%% supply the integrand: DW(\nabla u_h)*\grad w_h)^2
function val = integrandProj(x,y,curElem,lvl,p)

gradBasis = p.statics.stressBasis;
grad_h = p.statics.grad_h;
parents4e = p.level(lvl).enum.parents4e;
parent = parents4e(curElem);

evalGrad = grad_h(x,y,parent,lvl-1,p);
evalGradBasis = gradBasis(x,y,curElem,lvl,p);

evalGrad = reshape(evalGrad',[2 1 length(x)]);
val = matMul(evalGradBasis,evalGrad);
