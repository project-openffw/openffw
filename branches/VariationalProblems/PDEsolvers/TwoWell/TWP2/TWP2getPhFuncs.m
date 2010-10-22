function p = TWP2getPhFuncs(p)
% author: David Guenther, Lena Noack 

p.statics.sigma_h = @getSigma_h;

%% supply sigma_h = DW(\grad u_h)
function sigma_h = getSigma_h(x,y,curElem,lvl,p)
F1 = p.problem.F1;
F2 = p.problem.F2;

CONV = p.params.CONV;
if strcmp(CONV,'c')
nonLinearRegDer1 = p.problem.nonLinearRegDerA; %W'**1(F)
nonLinearRegDer2 = p.problem.nonLinearRegDerB; %W'**2(F)
else
nonLinearExactDer1 = p.problem.nonLinearExactDerA; %W'1(F)
nonLinearExactDer2 = p.problem.nonLinearExactDerB; %W'2(F)
end

grad_h = p.statics.grad_h;

curGrad = grad_h(x,y,curElem,lvl,p)';

if strcmp(CONV,'c')
%% sigma_h = DW**(F) = W'**1(F)F + W'**2(F)(F2,F)F2
sigma_h = curGrad*nonLinearRegDer1(curGrad(1,:)',curGrad(2,:)',curElem,lvl,p) +...
    F2*((curGrad(1,:)*F2(1)+curGrad(2,:)*F2(2))*nonLinearRegDer2(curGrad(1,:)',curGrad(2,:)',curElem,lvl,p));
else
%% sigma_h = DW(F) = W'1(F)(F-F1) + W'2(F)(F-F2)
% noch ändern, auf P2 anpassen
sigma_h = nonLinearExactDer1(curGrad(1),curGrad(2),curElem,lvl,p).*(curGrad-F1') +...
    nonLinearExactDer2(curGrad(1),curGrad(2),curElem,lvl,p).*(curGrad-F2');
end
sigma_h = ones(length(x),1) * sigma_h';
