function p = ODRTGOALgetJacobian(p)
% author: David Guenther
%% INPUT 
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
n4e = p.level(end).geom.n4e;

length4ed = p.level(end).enum.length4ed;
ed4e = p.level(end).enum.ed4e;
e4ed = p.level(end).enum.e4ed;
nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;

lvl = size(p.level,2);

degree = loadField('p.params','nonLinearExactIntegrateDegree',p,5);

%% compute the Jacobian DE(x)
BT = zeros(3,3,nrElems);
CT = zeros(1,3,nrElems);
DT = zeros(3,3,nrElems);

for curElem = 1:nrElems
	curEdges = ed4e(curElem,:);
	curEdgeLengths = length4ed(curEdges);
    
    signum = ones(1,3);
	I = find(e4ed(curEdges,2) == curElem);
	signum(I) = -1;    
    %divergence of basis functions on T,e.g. div(q) = sign*|E|/|T|
    div_qh = signum.*curEdgeLengths';
    
    D2W = integrate(n4e(curElem,:),lvl,degree,@integrand1,p);
%     D2W = intNonLinear(@getDer4Functional,sigma_h,stressBasis,stressBasis,[],...
%                             conjNonLinearFuncDer,conjNonLinearFuncSecDer,curElem,lvl,p);
    D3W = integrate(n4e(curElem,:),lvl,degree,@integrand2,p);
%     D3W = intNonLinear(@getSecDer4Functional,sigma_h,lambda1,stressBasis,stressBasis,...
%                             conjNonLinearFuncDer,conjNonLinearFuncSecDer,curElem,lvl,p);
    %calc int_T  D^2W(p_h)*q_h*r_h
    BT(:,:,curElem) = D2W;
   	%calc int_T v_h*div(q_h)
    CT(:,:,curElem) = div_qh;	
    
    %calc int_T  D^3W(p_h)*q_h*\lambda*r_h
    DT(:,:,curElem) = D3W;
end

[I,J] = localDoFtoGlobalDoF(ed4e,ed4e);
B = sparse(I,J,BT(:));
D = sparse(I,J,DT(:));

[I,J] = localDoFtoGlobalDoF(ed4e,(1:nrElems)');
C = sparse(I,J,CT(:),nrEdges,nrElems);

dummyZeros1 = sparse(nrEdges,nrElems);
dummyZeros2 = sparse(nrElems,nrElems);
dummyZeros3 = sparse(nrEdges,nrEdges);

jacobi = [  B-D,              dummyZeros1,     -B,         -C;
            dummyZeros1',  dummyZeros2,    C',      dummyZeros2;
            B          , C            ,  dummyZeros3,  dummyZeros1;
            C',  dummyZeros2, dummyZeros1', dummyZeros2];

%% OUTPUT 
p.level(end).jacobi = jacobi;

%% supply integrand
function val = integrand1(points,curElem,lvl,p)
% W''(|X|)/|X|^2*(X*Y)(X*Z) + W'(|X|)/|X|^3*( Y*Z*|X|^2 - (X*Y)(X*Z) )

DW = p.problem.conjNonLinearFuncDer;
D2W = p.problem.conjNonLinearFuncSecDer;
sigma_h = p.statics.sigma_h;
stressBasis = p.statics.stressBasis;

evalSigma = sigma_h(points,curElem,lvl,p);
evalBasis = stressBasis(points,curElem,lvl,p);

absSigma = ( evalSigma(:,1).^2 + evalSigma(:,2).^2 ).^(1/2);
evalDW = DW(absSigma,curElem,lvl,p);
evalD2W = D2W(absSigma,curElem,lvl,p);

evalSigma = reshape(evalSigma',[1 2 length(points(:,1))]);
YZ = matMul(evalBasis,permute(evalBasis,[2 1 3]));
XY = matMul(evalSigma,permute(evalBasis,[2 1 3]));
XZ = matMul(evalSigma,permute(evalBasis,[2 1 3]));
XYXZ = matMul(permute(XY,[2 1 3]),XZ);

if norm(absSigma) > 0
    term1 = matMul(reshape(evalD2W./absSigma.^2,[1 1 length(points(:,1))]),XYXZ);
    term2 = -matMul(reshape(evalDW./absSigma.^2,[1 1 length(points(:,1))]),XYXZ);
else
    term1 = 0;
    term2 = 0;
end

term3 = matMul(reshape(evalDW,[1 1 length(points(:,1))]),YZ);

val = term1 + term2 + term3;

%% supply integrand
function val = integrand2(points,curElem,lvl,p)
% computation of D^3W(x) = W'''(|x|)/|x|^3*(x*y)(x*z)(x*w) -...
% 2*W''(|x|)/|x|^3*(x*w)(x*z)(x*y) + W''(|x|)/|x|^2*( (x*y)(w*z)+(x*z)(w*y) )
% + W''(|x|)/|x|^2*(x*w)(y*z) - W''(|x|)/|x|^4*(x*z)(x*y)(x*w) - ...
% 3*W'(|x|)/|x|^2*(x*w)(y*z) + 3*W'(|x|)/|x|^4*(x*w)(x*z)(x*y) + ...
% W'(|x|)/|x|^3*(2*(y*z)(x*w)-(w*z)(x*y)-(x*z)(w*y))
%
% note that W depends on the method and we assume W''' = 0!!!!!!!!!!!!!


DW = p.problem.conjNonLinearFuncDer;
D2W = p.problem.conjNonLinearFuncSecDer;
sigma_h = p.statics.sigma_h;
lambda1 = p.statics.lambda1;
stressBasis = p.statics.stressBasis;

evalSigma = sigma_h(points,curElem,lvl,p);
evalLambda = lambda1(points,curElem,lvl,p);
evalBasis = stressBasis(points,curElem,lvl,p);

absSigma = ( evalSigma(:,1).^2 + evalSigma(:,2).^2 ).^(1/2);
evalDW = DW(absSigma,curElem,lvl,p);
evalD2W = D2W(absSigma,curElem,lvl,p);

evalSigma = reshape(evalSigma',[1 2 length(points(:,1))]);
evalLambda = reshape(evalLambda',[1 2 length(points(:,1))]);
YZ = matMul(evalLambda,permute(evalBasis,[2 1 3]));
XY = matMul(evalSigma,permute(evalLambda,[2 1 3]));
XZ = matMul(evalSigma,permute(evalBasis,[2 1 3]));
XYXZ = matMul(permute(XY,[2 1 3]),XZ);

YW = matMul(evalLambda,permute(evalBasis,[2 1 3]));
ZW = matMul(evalBasis,permute(evalBasis,[2 1 3]));
XY = matMul(evalSigma,permute(evalLambda,[2 1 3]));
XZ = matMul(evalSigma,permute(evalBasis,[2 1 3]));
XW = matMul(evalSigma,permute(evalBasis,[2 1 3]));

% evaluation of ansatz-function, dim = length(points(:,1))\times 2
% evalFunc = func(points,curElem,lvl,p);
% evaluation of the test-functions, 
% dim = squeeze([nrTestFunctions,2,length(points(:,1))])
% evalTestFunc1 = testFunc1(points,curElem,lvl,p);
% evalTestFunc2 = testFunc2(points,curElem,lvl,p);
% evalTestFunc3 = testFunc3(points,curElem,lvl,p);

% absFunc = (evalFunc(:,1).^2 + evalFunc(:,2).^2).^(1/2);

% evaluaton of the first nonLinear functions. dim = length(points(:,1))
% Note that nonLinear1 = nonLinear1(x)/x.
% evalNonLinear1 = nonLinear1(absFunc,curElem,lvl,p);
% Note that nonLinear2 = nonLinear2(x).
% evalNonLinear2 = nonLinear2(absFunc,curElem,lvl,p);

% if just one test function is supplied, we need to transform it in a 3-dim
% array
% sizeTestFunc1 = size(evalTestFunc1);
% if length(sizeTestFunc1) == 2
%     evalTestFunc1 = reshape(evalTestFunc1',[1 2 length(points(:,1))]);
% end
% sizeTestFunc2 = size(evalTestFunc2);
% if length(sizeTestFunc2) == 2
%     evalTestFunc2 = reshape(evalTestFunc2',[1 2 length(points(:,1))]);
% end
% sizeTestFunc3 = size(evalTestFunc3);
% if length(sizeTestFunc3) == 2
%     evalTestFunc3 = reshape(evalTestFunc3',[1 2 length(points(:,1))]);
% end

% evalFunc = reshape(evalFunc',[1 2 length(points(:,1))]);

% YZ = matMul(evalTestFunc1,permute(evalTestFunc2,[2 1 3]));
% YW = matMul(evalTestFunc1,permute(evalTestFunc3,[2 1 3]));
% ZW = matMul(evalTestFunc2,permute(evalTestFunc3,[2 1 3]));
% XY = matMul(evalFunc,permute(evalTestFunc1,[2 1 3]));
% XZ = matMul(evalFunc,permute(evalTestFunc2,[2 1 3]));
% XW = matMul(evalFunc,permute(evalTestFunc3,[2 1 3]));

XYZW = matMul(permute(XY,[2 1 3]),ZW);
XZYW = matMul(permute(XZ,[2 1 3]),YW);
XWYZ = matMul(permute(YZ,[2 1 3]),XW);
XZXY = matMul(permute(XZ,[2 1 3]),XY);
XWXZ = matMul(permute(XW,[2 1 3]),XZ);

XWXZXY = matMul(permute(XWXZ,[2 1 3]),XY);

if norm(absSigma) > 0
    % -2*W''(|x|)/|x|^3*(x*w)(x*z)(x*y)
    term1 = matMul(reshape(-2*evalD2W./absSigma.^3,[1 1 length(points(:,1))]),XWXZXY);
    % W''(|x|)/|x|^2*( (x*y)(w*z)+(x*z)(w*y) )
    term2 = matMul(reshape(evalD2W./absSigma.^2,[1 1 length(points(:,1))]),XYZW+XZYW);
    % W''(|x|)/|x|^2*(x*w)(y*z)
    term3 = matMul(reshape(evalD2W./absSigma.^2,[1 1 length(points(:,1))]),XWYZ);
    % - W''(|x|)/|x|^4*(x*z)(x*y)(x*w)
    term4 = matMul(reshape(-evalD2W./absSigma.^4,[1 1 length(points(:,1))]),XWXZXY);
    % -3*W'(|x|)/|x|^2*(x*w)(y*z)
    term5 = matMul(reshape(-3*evalDW./absSigma,[1 1 length(points(:,1))]),XWYZ);
    % 3*W'(|x|)/|x|^4*(x*w)(x*z)(x*y)
    term6 = matMul(reshape(3*evalDW./absSigma.^3,[1 1 length(points(:,1))]),XWXZXY);
    % W'(|x|)/|x|^3*(2*(y*z)(x*w)-(w*z)(x*y)-(x*z)(w*y))
    term7 = matMul(reshape(evalDW./absSigma.^2,[1 1 length(points(:,1))]),2*XWYZ-XYZW-XZYW);

    val = term1 + term2 + term3 + term4 + term5 + term6 + term7;
else
    val = zeros(3,3,length(points(:,1)));
end
