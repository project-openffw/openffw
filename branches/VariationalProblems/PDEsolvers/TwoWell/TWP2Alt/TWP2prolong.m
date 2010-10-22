function p = TWP2prolong(p,lvl)
% author: David Guenther, Lena Noack

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

%% supply integrand: DW(\nabla u_h)*\nabla w_h + 2(u_h - f)w_h
function val = integrandProj(x,y,curElem,lvl,p)

u_h = p.statics.u_h;
f = p.problem.f0;
basisU = p.statics.basisU;
evalUh = u_h(x,y,curElem,lvl,p);
evalF = f(x,y,curElem,lvl,p);
evalBasisU = basisU(x,y,curElem,lvl,p);

evalDif = reshape(2*(evalUh-evalF)',[1 1 length(x)]);
evalBasisU = reshape(evalBasisU,[3 1 length(x)]);

diffw = matMul(evalDif,evalBasisU);
rhs = diffw; %rhs = 2(u_h-f)w_h

sigma_h = p.statics.sigma_h;
stressBasis = p.statics.stressBasis;

evalSigma = sigma_h(x,y,curElem,lvl,p);
evalBasis = stressBasis(x,y,curElem,lvl,p);

evalSigma = reshape(evalSigma',[2 1 length(x)]);
val = matMul(evalBasis,evalSigma) + rhs; 
