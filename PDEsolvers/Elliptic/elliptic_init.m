function p = Elliptic_init(p)

% Set up integrands
p.statics.H1semiError = @H1semiErrorElliptic;
p.statics.L2error = @L2errorElliptic;


%%
function val = L2errorElliptic(x,y,curElem,lvl,p)

u_h = p.statics.u_h;
u_exact = p.problem.u_exact;

exactU =  u_exact(x, y,p)';
uh =  u_h(x, y,curElem,lvl,p);

val(1,1,:) = (uh - exactU).*(uh - exactU);


%%
function val = H1semiErrorElliptic(x,y,curElem,lvl,p)

sigma_h = p.statics.sigma_h;
gradU_exact = p.problem.gradU_exact;
kappa = p.problem.kappa;
% sigma_exact = kappa * gradU_exact;

curKappa = kappa(x,y,p);
if(ndims(curKappa) < 3)
    nrPoints = length(x)/length(curElem);
    curKappa = reshape(curKappa,nrPoints*length(curElem),4)';
    curKappa = reshape(curKappa,2,2,nrPoints,length(curElem));
end

curGradU_exact =  gradU_exact(x,y,p);
curGradU_exact = reshape (curGradU_exact',2,1,[]);
exactSigma = matMul(curKappa,curGradU_exact);
exactSigma = reshape(exactSigma,2,[])';

sigmah =  sigma_h(x,y,curElem,curKappa,lvl,p);

val(1,1,:) = sum( (sigmah - exactSigma).*(sigmah - exactSigma),2 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
