function p = ODRTestimate_oldP1Proj(p)
% author: David Guenther 

%% INPUT
u_D = p.problem.u_D;

c4n = p.level(end).geom.c4n;
n4e = p.level(end).geom.n4e;
Db = p.level(end).geom.Db;
grad4P1 = p.level(end).enum.grad4e;
area4e = p.level(end).enum.area4e;

nrElems = p.level(end).nrElems;
nrNodes = p.level(end).nrNodes;

lvl = size(p.level,2);

degree = loadField('p.params','nonLinearExactIntegrateDegree',p,19);

%% compute (\grad v_h,\grad w_h)-projection of DW*(p_h)
ST = zeros(3,3,nrElems);
bT = zeros(3,1,nrElems);

for curElem = 1:nrElems
    gradP1 = grad4P1(:,:,curElem);
    ST(:,:,curElem) = area4e(curElem)*gradP1*gradP1';
    DW = integrate(n4e(curElem,:),lvl,degree,@integrandProj,p);
    bT(:,:,curElem) = DW;
end

[I,J] = localDoFtoGlobalDoF(n4e,n4e);
S = sparse(I(:),J(:),ST(:));
dummy = n4e';
b = accumarray(dummy(:),bT(:));

DbNodes = unique(Db);

v = zeros(nrNodes,1);
v(DbNodes) = -u_D(c4n(DbNodes,1),c4n(DbNodes,2),p);
b = b - S*v;
freeNodes = setdiff(1:nrNodes,DbNodes);
v(freeNodes) = S(freeNodes,freeNodes)\b(freeNodes);

v4e = reshape(v(n4e'),[1 3 nrElems]);
gradVh = squeeze(matMul(v4e,grad4P1))';
p.level(end).P1proj = v;
p.level(end).gradP1proj4e = gradVh;

%% estimate the error by computing \eta=||DW*(p_h)-\grad v_h||_L^2

eta4T = zeros(nrElems,1);
for curElem = 1:nrElems    
    eta4T(curElem) = integrate(n4e(curElem,:),lvl,degree,@integrandError,p);
end

%% compose the error terms
estimatedError = sqrt( sum(eta4T) );

%% OUTPUT
p.level(end).etaT = sqrt(eta4T);
p.level(end).estimatedError = estimatedError;

%% supply the integrand: (DW*(p_h)-\grad v_h)^2
function val = integrandError(points,curElem,lvl,p)

gradVh = p.level(end).gradP1proj4e;
grad_h = p.statics.grad_h;

evalGrad = grad_h(points,curElem,lvl,p);

gradVh = repmat(gradVh(curElem,:),[length(points(:,1)),1]);
val = sum( (evalGrad-gradVh).^2,2 );
val = reshape(val,[1 1 length(points(:,1))]);

%% supply the integrand: DW*(p_h)*\grad w_h)^2
function val = integrandProj(points,curElem,lvl,p)

grad4P1 = p.level(end).enum.grad4e;
grad_h = p.statics.grad_h;

evalGrad = grad_h(points,curElem,lvl,p);
evalGradP1 = squeeze(grad4P1(:,:,curElem));

val = evalGradP1*evalGrad';
val = reshape(val,[3 1 length(points(:,1))]);