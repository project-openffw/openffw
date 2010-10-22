function p = ODP1DWRgetFuncValIM(p)
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
F_qh = integrate(n4e,lvl,degree,@integrandF,p);

I = dofU4e;
S = accumarray(I(:),ST(:));
rhs = accumarray(I(:),F_qh(:));

funcVal = rhs - S;

%% OUTPUT 
p.level(end).funcVal = funcVal(freeNodes);

%% supply integrand: D2W(\nabla u_h)*\nabla w_h*\nabla z_h
function val = integrand(x,y,curElem,lvl,p)
% W''(|X|)/|X|^2*(X*Y)(X*Z) + W'(|X|)/|X|^3*( Y*Z*|X|^2 - (X*Y)(X*Z) )

DW = p.problem.nonLinearExactDer;
D2W = p.problem.nonLinearExactSecDer;
grad_h = p.statics.grad_h;
DWRgrad_h = p.statics.DWRgrad_h;
stressBasis = p.statics.stressBasis;

evalGrad = grad_h(x,y,curElem,lvl,p);
evalDWRGrad = DWRgrad_h(x,y,curElem,lvl,p);
evalBasis = stressBasis(x,y,curElem,lvl,p);

absGrad = ( evalGrad(:,1).^2 + evalGrad(:,2).^2 ).^(1/2);
evalDW = DW(absGrad,curElem,lvl,p);
evalD2W = D2W(absGrad,curElem,lvl,p);

evalGrad = reshape(evalGrad',[1 2 length(x)]);
evalDWRGrad = reshape(evalDWRGrad',[1 2 length(x)]);
YZ = matMul(evalDWRGrad,permute(evalBasis,[2 1 3]));
XY = matMul(evalGrad,permute(evalDWRGrad,[2 1 3]));
XZ = matMul(evalGrad,permute(evalBasis,[2 1 3]));
XYXZ = matMul(permute(XY,[2 1 3]),XZ);

if norm(absGrad) > 0
    term1 = matMul(reshape(evalD2W./absGrad.^2,[1 1 length(x)]),XYXZ);
    term2 = -matMul(reshape(evalDW./absGrad.^2,[1 1 length(x)]),XYXZ);
else
    term1 = 0;
    term2 = 0;
end

term3 = matMul(reshape(evalDW,[1 1 length(x)]),YZ); %calculates elementwise matrix product

val = term1 + term2 + term3;








%% supply f_DWR * q_h (first goal function)
function val = integrandF(x,y,curElem,lvl,p)

%n4e = p.level(end).geom.n4e;
%c4n = p.level(end).geom.c4n;
%curNodes = n4e(curElem,:);
%Coords = c4n(curNodes,:);
stressBasis = p.statics.stressBasis; 
evalBasis = stressBasis(x,y,curElem,lvl,p);
%evalBasis = stressBasis(Coords(:,1),Coords(:,2),curElem,lvl,p);

f4e = p.level(end).DWRf4e;
evalF4e = f4e(curElem);
F4e = evalF4e*ones(2,length(x));
F4e = reshape(F4e,[1 2 length(x)]);

FQ = matMul(F4e,permute(evalBasis,[2 1 3]));

val = FQ;


%% supply f_DWR * v_h (first goal function)
function val = integrandFu(x,y,curElem,lvl,p)

Basis = p.statics.basisU; 
evalBasis = Basis(x,y,curElem,lvl,p);

f4e = p.level(end).DWRf4e;
evalF4e = f4e(curElem);
%F4e = evalF4e*ones(1,length(x));
%F4e = reshape(F4e,[1 1 length(x)]);

evalBasis=reshape(evalBasis,[1 3 length(x)]);
%FQ = matMul(F4e,permute(evalBasis,[2 1 3]));
val = evalF4e*evalBasis;


