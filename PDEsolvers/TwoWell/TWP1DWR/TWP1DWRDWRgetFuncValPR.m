%still to do: Goal function
function p = P1DWRgetFuncValPR(p)
%author: David Guenther, Lena Noack

%% INPUT
% load enumerated data
n4e = p.level(end).geom.n4e;
dofU4e = p.level(end).enum.dofU4e;
f4e = p.level(end).f4e;
lvl = size(p.level,2);
freeNodes = p.level(end).enum.freeNodes;
degree = loadField('p.params','nonLinearExactIntegrateDegree',p,5);

%% compute dual function value R(x)
ST = integrate(n4e,lvl,degree,@integrand,p);

Term_qh = integrate(n4e,lvl,degree,@evalDW,p);

I = dofU4e;
S = accumarray(I(:),ST(:));
rhs = accumarray(I(:),Term_qh(:));

funcVal = rhs - S;

%% OUTPUT 
p.level(end).funcVal = funcVal(freeNodes);

%% supply integrand: D2W(\nabla u_h)*\nabla w_h*\nabla z_h + 2 z_h w_h
function val = integrand(x,y,curElem,lvl,p)
F1 = p.problem.F1;
F2 = p.problem.F2;
CONV = p.params.CONV;

if strcmp(CONV,'c')
%% D2W**(F)[G,H] = W''**1(F)(G,H) + W''**2(F)(F,G)(F,H) + W''**3(F)(F2,H)(F2,G)
D2W1 = p.problem.nonLinearRegSecDerA;
D2W2 = p.problem.nonLinearRegSecDerB;
D2W3 = p.problem.nonLinearRegSecDerC;
else
% D2W(F)[G,H] = W''1(F)(G,H) + W''2(F)((F-F1,G)*(F-F2,H) + (F-F1,H)*(F-F2,G))
D2W1 = p.problem.nonLinearExactSecDerA;
D2W2 = p.problem.nonLinearExactSecDerB;
end


grad_h = p.statics.grad_h;
DWRu_h = p.statics.DWRu_h;
DWRgrad_h = p.statics.DWRgrad_h;
stressBasis = p.statics.stressBasis;
basisU = p.statics.basisU;

evalGrad = grad_h(x,y,curElem,lvl,p);
evalDWRu_h = DWRu_h(x,y,curElem,lvl,p);
evalDWRGrad = DWRgrad_h(x,y,curElem,lvl,p);
evalBasis = stressBasis(x,y,curElem,lvl,p);
evalBasisU = basisU(x,y,curElem,lvl,p);

grad1 = evalGrad(:,1);
grad2 = evalGrad(:,2);

evalD2W1 = D2W1(grad1,grad2,curElem,lvl,p);
evalD2W2 = D2W2(grad1,grad2,curElem,lvl,p);
if strcmp(CONV,'c')
evalD2W3 = D2W3(grad1,grad2,curElem,lvl,p);
end

evalGrad = reshape(evalGrad',[1 2 length(x)]);
evalDWRGrad = reshape(evalDWRGrad',[1 2 length(x)]);

GH = matMul(evalBasis,permute(evalDWRGrad,[2 1 3]));
evalGrad1 = reshape([grad1-F1(1) grad2-F1(2)]',[1 2 length(x)]);
evalGrad2 = reshape([grad1-F2(1) grad2-F2(2)]',[1 2 length(x)]);
F2vec = reshape([F2(1)*ones(size(grad1,1),1) F2(2)*ones(size(grad1,1),1)]',[1 2 length(x)]);
F1G = matMul(evalGrad1,permute(evalBasis,[2 1 3]));
F2G = matMul(evalGrad2,permute(evalBasis,[2 1 3]));
F1H = matMul(evalGrad1,permute(evalDWRGrad,[2 1 3]));
F2H = matMul(evalGrad2,permute(evalDWRGrad,[2 1 3]));
FG = matMul(evalGrad,permute(evalBasis,[2 1 3]));
FH = matMul(evalGrad,permute(evalDWRGrad,[2 1 3]));
GF2 = matMul(F2vec,permute(evalBasis,[2 1 3]));
HF2 = matMul(F2vec,permute(evalDWRGrad,[2 1 3]));
F1GF2H = matMul(permute(F1G,[2 1 3]),F2H);
F1HF2G = matMul(permute(F1H,[2 1 3]),F2G);
FGFH = matMul(permute(FG,[2 1 3]),FH);
GF2HF2 = matMul(permute(GF2,[2 1 3]),HF2);

%% D2W**(F)[G,H] = W''**1(F)(G,H) + W''**2(F)(F,G)(F,H) + W''**3(F)(F2,G)(F2,H)
% D2W(F)[G,H] = W''1(F)(G,H) + W''2(F)((F-F1,G)*(F-F2,H) + (F-F1,H)*(F-F2,G))

if strcmp(CONV,'c')
term1 = matMul(reshape(evalD2W1,[1 1 length(x)]),GH)+matMul(reshape(evalD2W2,[1 1 length(x)]),FGFH);
term2 = matMul(reshape(evalD2W3,[1 1 length(x)]),GF2HF2);
else
term1 = matMul(reshape(evalD2W1,[1 1 length(x)]),GH);
term2 = matMul(reshape(evalD2W2,[1 1 length(x)]),F1GF2H+F1HF2G);
end

val = term1 + term2 + 2*matMul(reshape(evalBasisU,[3 1 length(x)]),reshape(evalDWRu_h,[1 1 length(x)]));


function val = evalDW(x,y,curElem,lvl,p)

sigma_h = p.statics.sigma_h;
DW_Vh = sigma_h(x,y,curElem,lvl,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F1 = p.problem.F1;
F2 = p.problem.F2;
CONV = p.params.CONV;

if strcmp(CONV,'c')
%% D2W**(F) = W''**1(F) + W''**2(F)(F,F) + W''**3(F)(F2,F2)
D2W1 = p.problem.nonLinearRegSecDerA;
D2W2 = p.problem.nonLinearRegSecDerB;
D2W3 = p.problem.nonLinearRegSecDerC;
else
% D2W(F) = W''1(F) + 2*W''2(F)(F-F1,F-F2)
D2W1 = p.problem.nonLinearExactSecDerA;
D2W2 = p.problem.nonLinearExactSecDerB;
end

grad_h = p.statics.grad_h;
evalGrad = grad_h(x,y,curElem,lvl,p);

grad1 = evalGrad(:,1);
grad2 = evalGrad(:,2);

evalGrad = reshape(evalGrad',[1 2 length(x)]);

evalD2W1 = D2W1(grad1,grad2,curElem,lvl,p);
evalD2W2 = D2W2(grad1,grad2,curElem,lvl,p);
if strcmp(CONV,'c')
evalD2W3 = D2W3(grad1,grad2,curElem,lvl,p);
end

evalGrad1 = reshape([grad1-F1(1) grad2-F1(2)]',[1 2 length(x)]);
evalGrad2 = reshape([grad1-F2(1) grad2-F2(2)]',[1 2 length(x)]);
F2vec = reshape([F2(1)*ones(size(grad1,1),1) F2(2)*ones(size(grad1,1),1)]',[1 2 length(x)]);
F1F2 = matMul(evalGrad1,permute(evalGrad2,[2 1 3])); % (F-F1,F-F2)
FF = matMul(evalGrad,permute(evalGrad,[2 1 3]));     % (F,F)
F2F2 = matMul(F2vec,permute(F2vec,[2 1 3]));         % (F2,F2)


if strcmp(CONV,'c')
%% D2W**(F) = W''**1(F) + W''**2(F)(F,F) + W''**3(F)(F2,F2)
term1 = reshape(evalD2W1,[1 1 length(x)]) + matMul(reshape(evalD2W2,[1 1 length(x)]),FF);
term2 = matMul(reshape(evalD2W3,[1 1 length(x)]),F2F2);
else
% D2W(F) = W''1(F) + 2*W''2(F)(F-F1,F-F2)
term1 = reshape(evalD2W1,[1 1 length(x)]);
term2 = 2*matMul(reshape(evalD2W2,[1 1 length(x)]),F1F2);
end

D2W_Vh = term1 + term2;                        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DW_Vh = reshape(sum(DW_Vh,2),[1 1 length(x)]);
val = 2*matMul(DW_Vh,permute(D2W_Vh,[2 1 3]));       
val = permute(val,[2 1 3]);

val = matMul(reshape(ones(3,length(x)),[3 1 length(x)]),val.^2);
