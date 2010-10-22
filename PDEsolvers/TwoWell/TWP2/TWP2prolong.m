function p=TWP2prolong(p,lvl)
% author: David Guenther, Lena Noack

%% INPUT
nrNodes = p.level(end).nrNodes;
nrEdges = p.level(end).nrEdges;
nrElems = p.level(lvl).nrElems;
parents4e = p.level(lvl).enum.parents4e;

%% prolongation
n4e = p.level(end).geom.n4e;
lvl = size(p.level,2);
degree = loadField('p.params','nonLinearExactIntegrateDegree',p,19);

x0 = zeros(nrNodes+nrEdges,1);
if lvl > 2
    q = initFFW('P2',p.params.problem.name,'uniform',p.level(end).nrDoF,'elliptic','direct','redGreenBlue','STUBestimate',false);
    q.params.errorOutput = false;
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

%% supply the integrand: DW(\nabla u_h)*\nabla w_h + 2(u_h - f)w_h
function val = integrandProj(x,y,curElem,lvl,p)

gradBasis = p.statics.stressBasis;
grad_h = p.statics.grad_h;
parents4e = p.level(lvl).enum.parents4e;
parent = parents4e(curElem);

u_h = p.statics.u_h;
f = p.problem.f0;
evalU = u_h(x,y,parent,lvl-1,p);
evalF = f(x,y,parent,lvl-1,p);

evalGrad = grad_h(x,y,parent,lvl-1,p);
evalGradBasis = gradBasis(x,y,curElem,lvl,p);

evalGrad = reshape(evalGrad',[2 1 length(x)]);

% second part 
basisU = p.statics.basisU;
evalBasisU = basisU(x,y,curElem,lvl,p);

evalBasisU=reshape(evalBasisU,[6 1 length(x)]);
evalBasisU2=permute(evalBasisU,[2 1 3]);
evalBasisU2=reshape(evalBasisU2,[1 6 length(x)]);

evalDif = reshape(2*(evalU-evalF)',[1 1 length(x)]);
rhs = matMul(evalDif,evalBasisU); %rhs = 2(u_h-f)w_h


val = matMul(evalGradBasis,evalGrad)+ rhs;