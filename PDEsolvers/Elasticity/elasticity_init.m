function p = Elasticity_init(p)

    p.statics.H1semiError = @H1semiErrorElasticity;
    p.statics.L2error = @L2errorElasticity;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function val = L2errorElasticity(x,y,curElem,lvl,p)

u_h = p.statics.u_h;
u_exact = p.problem.u_exact;

exactU =  u_exact(x,y,p)';

uh =  u_h(x,y,curElem,lvl,p);
val(1,1,:) = sum((uh - exactU).*(uh - exactU),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = H1semiErrorElasticity(x,y,curElem,lvl,p)

sigma_h = p.statics.sigma_h;
sigma_exact = p.problem.sigma_exact;

exactSigma =  sigma_exact(x,y,p);
sigmah =  sigma_h(x,y,curElem,lvl,p);

val(1,1,:) = sum(sum( (sigmah - exactSigma).*(sigmah - exactSigma),2 ),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
