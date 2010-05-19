function p = optimalDesign_init(p)
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

p.statics.volumeFraction = @getVolumeFraction;
p.statics.energy_h = @getEnergy_h;

p.statics.L2errorDisplacement = @L2errorDisplacement;
p.statics.L2errorGradU = @L2errorGradU;
p.statics.L2errorPhminusP0 = @L2errorPhminusP0;
p.statics.L2errorPminusNonLinear = @L2errorPminusNonLinear;
p.statics.errorEnergy = @errorEnergy;

%% supply discrete energy functional
function val = getEnergy_h(points,curElem,lvl,p)

grad_h = p.statics.grad_h;
nonLinearExact = p.problem.nonLinearExact;
u_h = p.statics.u_h;
f = p.problem.f;
area4e = p.level(lvl).enum.area4e;
lambda = p.problem.lambda;
mu1 = p.problem.mu1;
mu2 = p.problem.mu2;
distributeParam = loadField('p.params','distributeParam',p,0.5);

gradU = grad_h(points,curElem,lvl,p);
absGrad = (gradU(:,1).^2+gradU(:,2).^2).^(1/2);
evalNonLinear = nonLinearExact(absGrad,curElem,lvl,p);

evalU = u_h(points,curElem,lvl,p);
evalF = f(points,curElem,lvl,p);
% evalF = f(x,y,p);

optDesign = lambda*(distributeParam*mu1 + (1-distributeParam)*mu2);
% optDesign = 0;
val = evalNonLinear - evalF.*evalU + optDesign;
val = reshape(val,[1 1 length(points(:,1))]);

%% supply volume fractions
function val = getVolumeFraction(points,curElem,lvl,p)

grad_h = p.statics.grad_h;

gradU = grad_h(points,curElem,lvl,p);
absGrad = (gradU(:,1).^2+gradU(:,2).^2).^(1/2);

t1 = p.problem.t1;
t2 = p.problem.t2;

val = zeros(length(points(:,1)),1);
% val(absGrad < t1) = (t1 - absGrad(absGrad < t1))/(t2-t1);
index = find(absGrad >= t1 & absGrad <= t2);
% val(index) = (t2-absGrad(index))/(t2-t1);
val(index) = (absGrad(index)-t1)/(t2-t1);
% val(absGrad > t2) = (absGrad(absGrad > t2) - t2)/(t2-t1);
val(absGrad>t2) = 1;

%% supply the error ||u-u_h||_L^2
function val = L2errorDisplacement(points,curElem,lvl,p)

u_h = p.statics.u_h;
u_exact = p.problem.u_exact;

evalUexact = u_exact(points,p);
evalUh = u_h(points,curElem,lvl,p);

val = (evalUh - evalUexact).*(evalUh - evalUexact);
val = reshape(val,[1 1 length(points(:,1))]);

%% supply the error ||\grad(u-u_h)||_L^2
function val = L2errorGradU(points,curElem,lvl,p)

gradU_h = p.statics.grad_h;
gradU_exact = p.problem.gradU_exact;

evalGradUexact = gradU_exact(points,p);
evalGradUh = gradU_h(points,curElem,lvl,p);

val = sum((evalGradUh - evalGradUexact).*(evalGradUh - evalGradUexact),2);
val = reshape(val,[1 1 length(points(:,1))]);

%% supply the error ||p_{h,\epsilon}-p_\epsilon||_L^2
function val = L2errorPminusNonLinear(points,curElem,lvl,p)

sigma_Eps = p.problem.sigmaEps;
sigma_h = p.statics.sigma_h;

sigmaEps =  sigma_Eps(points,curElem,lvl,p);
sigmah =  sigma_h(points,curElem,lvl,p);

val = sum( (sigmah - sigmaEps).*(sigmah - sigmaEps),2 );
val = reshape(val,[1 1 length(points(:,1))]);

%% supply the error ||p_{h,\epsilon)-p||_L^2 with p = DW(\grad u)
function val = L2errorPhminusP0(points,curElem,lvl,p)

sigma_0 = p.problem.sigma0;
sigma_h = p.statics.sigma_h;

sigma0 = sigma_0(points,curElem,lvl,p);
sigmah = sigma_h(points,curElem,lvl,p);

val = sum( (sigmah - sigma0).*(sigmah - sigma0),2 );
val = reshape(val,[1 1 length(points(:,1))]);

%% supply the error in the energy-functional \int_T |W(\grad u_h)-W(\grad u)-f(u_h-u)| 
function val = errorEnergy(points,curElem,lvl,p)

u_exact = p.problem.u_exact;
% I2u_h = p.statics.I2u_h;
gradU_exact = p.problem.gradU_exact;
% I2p_h = p.statics.I2p_h;
nonLinearExact = p.problem.nonLinearExact;
conjNonLinearFuncDer = p.problem.conjNonLinearFuncDer;
f = p.problem.f;

u_h = p.statics.u_h;
grad_h = p.statics.grad_h;

evalUexact = u_exact(points,p);
evalGradUexact = gradU_exact(points,p);
% evalUexact = I2u_h(x,y,curElem,lvl,p);
% evalI2p_h = I2p_h(x,y,curElem,lvl,p);

% absI2p_h = ( evalI2p_h(:,1).^2 + evalI2p_h(:,2).^2 ).^(1/2);
% evalConjNonLinear = conjNonLinearFuncDer(absI2p_h,curElem,lvl,p).*absI2p_h;

% evalGradUexact = (evalConjNonLinear*[1,1]).*evalI2p_h;

absGradExact = ( evalGradUexact(:,1).^2 + evalGradUexact(:,2).^2 ).^(1/2);

evalUh = u_h(points,curElem,lvl,p);
evalGradh = grad_h(points,curElem,lvl,p);

absGradh = ( evalGradh(:,1).^2 + evalGradh(:,2).^2 ).^(1/2);

evalNonLinearExact = nonLinearExact(absGradExact,curElem,lvl,p);
evalNonLinearDiscrete = nonLinearExact(absGradh,curElem,lvl,p);

% evalF = f(x,y,p);
evalF = f(points,curElem,lvl,p);

val = abs((evalNonLinearDiscrete - evalNonLinearExact) - ...
        evalF.*(evalUh-evalUexact));

val = reshape(val,[1 1 length(points(:,1))]);
