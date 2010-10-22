function p = TWP1estimate_JumpCC(p)
%author: Lena Noack
% ( h_T^2 ||f +div sigma_h||^2_(L^4/3) + 0.5 sum h_E ||[sigma_h]nu||_{L^4/3} )^{1/2}

%% INPUT
length4ed = p.level(end).enum.length4ed;
n4e = p.level(end).geom.n4e;
n4ed = p.level(end).enum.n4ed;
ed4e = p.level(end).enum.ed4e;
curLvl = length(p.level);
degree = loadField('p.params','rhsIntegrateExactDegree',p,5);

grad_h = p.statics.grad_h;  % grad_h=p_h
e4ed = p.level(end).enum.e4ed;
length4ed = p.level(end).enum.length4ed;
midPoint4e = p.level(end).enum.midPoint4e;
nrElems = p.level(end).nrElems;
lvl = length(p.level);

normals4ed = p.level(lvl).enum.normals4ed;
nrEdges = p.level(end).nrEdges;

h_T = max(length4ed(ed4e),[],2);

mu_T = integrate(n4e,curLvl,2*degree,@Residuum,p);
mu_E = integrate(n4ed,curLvl,degree,@integrand,p);
mu = sqrt( h_T.^2.*mu_T + 1/2*h_T.*sum(mu_E(ed4e),2) );

%% OUTPUT 
p.level(end).etaT = mu;
p.level(end).estimatedError = (norm(mu,2)).^(3/4);
%p.level(end).estimatedError = norm(mu,2);

function val = Residuum(x,y,curElem,curLvl,p)
f = p.problem.f;
curf     = f(x,y,curElem,curLvl,p);

F1 = p.problem.F1;
F2 = p.problem.F2;
CONV = p.params.CONV;

if strcmp(CONV,'c')
%% D2W**(F)[G,H] = W''**1(F) + W''**2(F)(F,F) + W''**3(F)(F2,F2)
D2W1 = p.problem.nonLinearRegSecDerA;
D2W2 = p.problem.nonLinearRegSecDerB;
D2W3 = p.problem.nonLinearRegSecDerC;
else
% D2W(F)[G,H] = W''1(F) + W''2(F)*2*(F-F1,F-F2)
D2W1 = p.problem.nonLinearExactSecDerA;
D2W2 = p.problem.nonLinearExactSecDerB;
end

grad_h = p.statics.grad_h;
stressBasis = p.statics.stressBasis;
basisU = p.statics.basisU;

evalGrad = grad_h(x,y,curElem,curLvl,p);

grad1 = evalGrad(:,1);
grad2 = evalGrad(:,2);
evalGrad = reshape(evalGrad',[1 2 length(x)]);

evalD2W1 = D2W1(grad1,grad2,curElem,curLvl,p);
evalD2W2 = D2W2(grad1,grad2,curElem,curLvl,p);
if strcmp(CONV,'c')
evalD2W3 = D2W3(grad1,grad2,curElem,curLvl,p);
end

evalGrad1 = reshape([grad1-F1(1) grad2-F1(2)]',[1 2 length(x)]);
evalGrad2 = reshape([grad1-F2(1) grad2-F2(2)]',[1 2 length(x)]);
F2vec = reshape([F2(1)*ones(size(grad1,1),1) F2(2)*ones(size(grad1,1),1)]',[1 2 length(x)]);

F1F2 = matMul(evalGrad1,permute(evalGrad2,[2 1 3]));
FF = matMul(evalGrad,permute(evalGrad,[2 1 3]));
F2F2 = matMul(F2vec,permute(F2vec,[2 1 3]));

if strcmp(CONV,'c')
term1 = reshape(evalD2W1,[1 1 length(x)])+matMul(reshape(evalD2W2,[1 1 length(x)]),FF);
term2 = matMul(reshape(evalD2W3,[1 1 length(x)]),F2F2);
else
term1 = reshape(evalD2W1,[1 1 length(x)]);
term2 = 2*matMul(reshape(evalD2W2,[1 1 length(x)]),F1F2);
end

evalD2W = term1 + term2;






residuum = - curf(:);

evalD2W=reshape(evalD2W,[length(x) 1]);
sumRes = residuum + evalD2W;

val(1,:,:) = (sumRes.^2)';

function val = integrand(x,y,curEdge,lvl,p)

sigma_h = p.statics.sigma_h; %sigma_h = DW(p_h)
e4ed = p.level(lvl).enum.e4ed;
normals4ed = p.level(lvl).enum.normals4ed;

elems = e4ed(curEdge,:);
normal = normals4ed(curEdge,:);

evalSigma1 = sigma_h(x,y,elems(1),lvl,p)*normal'; %sigma_h|T_1 * nu_E

if elems(2) ~= 0
    % E is an interior edge
    evalSigma2 = sigma_h(x,y,elems(2),lvl,p)*normal'; %sigma_h|T_2 * nu_E
else
    % E is a Dirichlet Edge
    evalSigma2 = evalSigma1;
end

val = zeros(1,1,length(x));
val(1,1,:) = sum((evalSigma1 - evalSigma2).^2,2); %sum( (sigma_h|T_1 - sigma_h|T_2 ) * nu_E )^2
