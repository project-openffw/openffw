function p = ODP2getPhFuncs(p)
% author: David Guenther 

p.statics.sigma_h = @getSigma_h;

%% supply p_h = DW(\grad u_h)
function sigma_h = getSigma_h(points,curElem,lvl,p)

nonLinearExactDer = p.problem.nonLinearExactDer;
grad_h = p.statics.grad_h;

evalGrad = grad_h(points,curElem,lvl,p);
absGrad = ( evalGrad(:,1).^2 + evalGrad(:,2).^2 ).^(1/2);

sigma_h = (nonLinearExactDer(absGrad,curElem,lvl,p)*[1 1]).*evalGrad;
