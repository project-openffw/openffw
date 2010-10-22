function p = AbsValue_init(p)

p.statics.volumeFraction = @getVolumeFraction; %rausschmeiﬂen
p.statics.energy_h = @getEnergy_h;

p.statics.L2errorDisplacement = @L2errorDisplacement;
p.statics.L2errorGradU = @L2errorGradU;
p.statics.L2errorPhminusP0 = @L2errorPhminusP0;
p.statics.L2errorPhminusP0Square = @L2errorPhminusP0Square;
p.statics.L2errorPminusNonLinear = @L2errorPminusNonLinear;
p.statics.errorEnergy = @errorEnergy;

%% supply discrete energy functional
function val = getEnergy_h(x,y,curElem,lvl,p)

grad_h = p.statics.grad_h;
nonLinearExact = p.problem.nonLinearExact;
u_h = p.statics.u_h;
f = p.problem.f;
area4e = p.level(lvl).enum.area4e;
pv = p.problem.pv;
distributeParam = loadField('p.params','distributeParam',p,0.5);

gradU = grad_h(x,y,curElem,lvl,p);
absGrad = (gradU(:,1).^2+gradU(:,2).^2).^(1/2);
evalNonLinear = nonLinearExact(absGrad,curElem,lvl,p);

evalU = u_h(x,y,curElem,lvl,p);
evalF = f(x,y,curElem,lvl,p);
% evalF = f(x,y,p);

val = evalNonLinear - evalF.*evalU;
val = reshape(val,[1 1 length(x)]);

%% supply volume fractions
function val = getVolumeFraction(x,y,curElem,lvl,p)

val = 1;

%% supply the error ||u-u_h||_L^2
function val = L2errorDisplacement(x,y,curElem,lvl,p)

u_h = p.statics.u_h;
u_exact = p.problem.u_exact;

evalUexact = u_exact(x,y,p);
evalUh = u_h(x,y,curElem,lvl,p);

val = (evalUh - evalUexact).*(evalUh - evalUexact);
val = reshape(val,[1 1 length(x)]);

%% supply the error ||\grad(u-u_h)||_L^2
function val = L2errorGradU(x,y,curElem,lvl,p)

gradU_h = p.statics.grad_h;
gradU_exact = p.problem.gradU_exact;

evalGradUexact = gradU_exact(x,y,p);
evalGradUh = gradU_h(x,y,curElem,lvl,p);

val = sum((evalGradUh - evalGradUexact).*(evalGradUh - evalGradUexact),2);
val = reshape(val,[1 1 length(x)]);

%% supply the error ||p_{h,\epsilon}-p_\epsilon||_L^2
function val = L2errorPminusNonLinear(x,y,curElem,lvl,p)

sigma_Eps = p.problem.sigmaEps;
sigma_h = p.statics.sigma_h;

sigmaEps =  sigma_Eps(x,y,curElem,lvl,p);
sigmah =  sigma_h(x,y,curElem,lvl,p);

val = sum( (sigmah - sigmaEps).*(sigmah - sigmaEps),2 );
val = reshape(val,[1 1 length(x)]);

%% supply the error ||p_h-p||_L^2 with p = DW(\grad u)
function val = L2errorPhminusP0(x,y,curElem,lvl,p)

sigma_0 = p.problem.sigma0;
sigma_h = p.statics.sigma_h;

sigma0 = sigma_0(x,y,curElem,lvl,p);
sigmah = sigma_h(x,y,curElem,lvl,p);

val = sum( (sigmah - sigma0).*(sigmah - sigma0),2 );
val = reshape(val,[1 1 length(x)]);

%% supply the error ||p_h-p||^2_L^2 with p = DW(\grad u)
function val = L2errorPhminusP0Square(x,y,curElem,lvl,p)

%fprintf('SquareInit \n')
p.params.square=1;
sigma_0 = p.problem.sigma0;
sigma_h = p.statics.sigma_h;

sigma0 = sigma_0(x,y,curElem,lvl,p);
sigmah = sigma_h(x,y,curElem,lvl,p);

val = sum( (sigmah - sigma0).*(sigmah - sigma0),2 );
val = reshape(val,[1 1 length(x)]);

%% supply the error in the energy-functional \int_T |W(\grad u_h)-W(\grad u)-f(u_h-u)| 
function val = errorEnergy(x,y,curElem,lvl,p)

u_exact = p.problem.u_exact;
% I2u_h = p.statics.I2u_h;
gradU_exact = p.problem.gradU_exact;
% I2p_h = p.statics.I2p_h;
nonLinearExact = p.problem.nonLinearExact;
nonLinearFuncDer = p.problem.nonLinearExactDer;
f = p.problem.f;

u_h = p.statics.u_h;
grad_h = p.statics.grad_h;

evalUexact = u_exact(x,y,p);
evalGradUexact = gradU_exact(x,y,p);

absGradExact = ( evalGradUexact(:,1).^2 + evalGradUexact(:,2).^2 ).^(1/2);

evalUh = u_h(x,y,curElem,lvl,p);
evalGradh = grad_h(x,y,curElem,lvl,p);

absGradh = ( evalGradh(:,1).^2 + evalGradh(:,2).^2 ).^(1/2);

evalNonLinearExact = nonLinearExact(absGradExact,curElem,lvl,p);
evalNonLinearDiscrete = nonLinearExact(absGradh,curElem,lvl,p);

% evalF = f(x,y,p);
evalF = f(x,y,curElem,lvl,p);

val = abs((evalNonLinearDiscrete - evalNonLinearExact) - ...
        evalF.*(evalUh-evalUexact));

val = reshape(val,[1 1 length(x)]);
