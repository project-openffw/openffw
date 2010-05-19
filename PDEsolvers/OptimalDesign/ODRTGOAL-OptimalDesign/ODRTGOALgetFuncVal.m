function p = ODRTGOALgetFuncVal(p)
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
n4ed = p.level(end).enum.n4ed;
DbEd = p.level(end).enum.DbEd;

nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;

x = p.level(end).x0;
lvl = size(p.level,2);

degree = loadField('p.params','nonLinearExactIntegrateDegree',p,5);

%% compute the function value E(x)

BT = zeros(1,3,nrElems);
CT = zeros(1,3,nrElems);
DT = zeros(1,1,nrElems);
ET = zeros(1,3,nrElems);
FT = zeros(1,3,nrElems);
GT = zeros(1,1,nrElems);

for curElem = 1:nrElems
	curEdges = ed4e(curElem,:);
	curEdgeLengths = length4ed(curEdges);
    curX4ed = x(curEdges);
    curX4e = x(nrEdges + curElem);
    curLambda1 = x(nrEdges + nrElems + curEdges);
    curLambda2 = x(nrEdges + nrElems + nrEdges + curElem);
    
    signum = ones(1,3);
	I = find(e4ed(curEdges,2) == curElem);
	signum(I) = -1;    
    %divergence of basis functions on T,e.g. div(q) = sign*|E|/|T|
    div_qh = signum.*curEdgeLengths';
    
    %calc int_T  DW(p_h)*q_h
    DW = integrate(n4e(curElem,:),lvl,degree,@integrand1,p);

    %calc int_T  D^2W(p_h)*\lambda1*q_h
    D2W = integrate(n4e(curElem,:),lvl,degree,@integrand2,p);
    D2W = D2W(:)';

    BT(:,:,curElem) = DW;
   	%calc int_T u_h*div(q_h)
    CT(:,:,curElem) = curX4e*div_qh;	
    %calc int_T v_h*div(p_h)
    DT(:,:,curElem) = curX4ed'*div_qh';	
    
    ET(:,:,curElem) = DW - D2W;
    FT(:,:,curElem) = curLambda2*div_qh;
    
    GT(:,:,curElem) = curLambda1'*div_qh';	
end

I = ed4e';
B = accumarray(I(:),BT(:));
C = accumarray(I(:),CT(:));
E = accumarray(I(:),ET(:));
F = accumarray(I(:),FT(:));

I = 1:nrElems;
D = accumarray(I(:),DT(:));
G = accumarray(I(:),GT(:));

boundary = zeros(nrEdges + nrElems,1);
boundary(DbEd) = integrate(n4ed(DbEd,:),lvl,degree,@integrandBoundary,p);

b = -boundary;
f4e = p.level(end).f4e;
b(nrEdges+1:end) = f4e;

% b = -boundary;
% f4e = p.level(end).f4e;
% b(nrEdges+1:end) = f4e;

funcVal1 = [B + C;
               D]   +   b;

funcVal2 = [E - F;
            G];

funcVal = [funcVal2;
           funcVal1];

%% OUTPUT 
p.level(end).funcVal = funcVal;

%% supply integrand DW^*(p_h)*q_h
function val = integrand1(points,curElem,lvl,p)
% W'(|X|)/|X|*X*Y

grad_h = p.statics.grad_h;
stressBasis = p.statics.stressBasis;

evalGrad = grad_h(points,curElem,lvl,p);
evalGrad = reshape(evalGrad',[2 1 length(points(:,1))]);

evalBasis = stressBasis(points,curElem,lvl,p);

val = matMul(evalBasis,evalGrad);

%% supply integrand
function val = integrand2(points,curElem,lvl,p)
% W''(|X|)/|X|^2*(X*Y)(X*Z) + W'(|X|)/|X|^3*( Y*Z*|X|^2 - (X*Y)(X*Z) )

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

if norm(absSigma) > 0
    term1 = matMul(reshape(evalD2W./absSigma.^2,[1 1 length(points(:,1))]),XYXZ);
    term2 = -matMul(reshape(evalDW./absSigma.^2,[1 1 length(points(:,1))]),XYXZ);
else
    term1 = 0;
    term2 = 0;
end

term3 = matMul(reshape(evalDW,[1 1 length(points(:,1))]),YZ);

val = term1 + term2 + term3;

%% supply integrand for Dirichlet boundary
function val = integrandBoundary(points,curEdge,lvl,p)

u_D = p.problem.u_D;
stressBasis = p.statics.stressBasis;
e4ed = p.level(lvl).enum.e4ed;
ed4e = p.level(lvl).enum.ed4e;
normals4ed = p.level(lvl).enum.normals4ed;

curElem = e4ed(curEdge,1);
edges = ed4e(curElem,:);
index = find(edges == curEdge);
normal = normals4ed(curEdge,:);

evalU_D = u_D(points,p);

evalStressBasis = stressBasis(points,curElem,lvl,p);
curBasisFunc = squeeze(evalStressBasis(index,:,:))';
basis4normal = curBasisFunc*normal';

val = evalU_D.*basis4normal;
val = reshape(val,[1 1 length(points(:,1))]);
