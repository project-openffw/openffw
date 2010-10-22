function p = TwoWell_init(p)

p.statics.volumeFraction = @getVolumeFraction;
p.statics.energy_h = @getEnergy_h;

p.statics.L43errorDisplacement = @L43errorDisplacement;
p.statics.L43errorGradU = @L43errorGradU;
p.statics.L43errorPhminusP0 = @L43errorPhminusP0;
p.statics.L43errorPhminusP0Square = @L43errorPhminusP0Square;
p.statics.L43errorPminusNonLinear = @L43errorPminusNonLinear;
p.statics.errorEnergy = @errorEnergy;

%% supply discrete energy functional
function val = getEnergy_h(x,y,curElem,lvl,p)

grad_h = p.statics.grad_h;
CONV = p.params.CONV;
if strcmp(CONV,'c')
nonLinearExact = p.problem.nonLinearReg;
else
nonLinearExact = p.problem.nonLinearExact;
end
u_h = p.statics.u_h;
f = p.problem.f0;
area4e = p.level(lvl).enum.area4e;
distributeParam = loadField('p.params','distributeParam',p,0.5);

gradU = grad_h(x,y,curElem,lvl,p);
evalNonLinear = nonLinearExact(gradU(:,1),gradU(:,2),curElem,lvl,p);

evalU = u_h(x,y,curElem,lvl,p);
evalF = f(x,y,curElem,lvl,p);
% evalF = f(x,y,p);

val = evalNonLinear + (evalF-evalU).*(evalF-evalU);
val = reshape(val,[1 1 length(x)]);

%% supply volume fractions
function val = getVolumeFraction(x,y,curElem,lvl,p)

val = 1;

%% supply the error ||u-u_h||_L^4/3
function val = L43errorDisplacement(x,y,curElem,lvl,p)

u_h = p.statics.u_h;
u_exact = p.problem.u_exact;

evalUexact = u_exact(x,y,p);
evalUh = u_h(x,y,curElem,lvl,p);

val = (evalUh - evalUexact).^(4/3);
val = reshape(val,[1 1 length(x)]);

%% supply the error ||\grad(u-u_h)||_L^4/3
function val = L43errorGradU(x,y,curElem,lvl,p)

gradU_h = p.statics.grad_h;
gradU_exact = p.problem.gradU_exact;

evalGradUexact = gradU_exact(x,y,p);
evalGradUh = gradU_h(x,y,curElem,lvl,p);

val = sum((evalGradUh - evalGradUexact).^(4/3),2);
val = reshape(val,[1 1 length(x)]);

%% supply the error ||p_{h,\epsilon}-p_\epsilon||_L^4/3
function val = L43errorPminusNonLinear(x,y,curElem,lvl,p)

sigma_Eps = p.problem.sigmaEps;
sigma_h = p.statics.sigma_h;

sigmaEps =  sigma_Eps(x,y,curElem,lvl,p);
sigmah =  sigma_h(x,y,curElem,lvl,p);

val = sum( (sigmah - sigmaEps).^(4/3),2 );
val = reshape(val,[1 1 length(x)]);

%% supply the error ||p_h-p||_L^(4/3) with p = DW(\grad u)
function val = L43errorPhminusP0(x,y,curElem,lvl,p)

sigma_0 = p.problem.sigma0;
sigma_h = p.statics.sigma_h;

sigma0 = sigma_0(x,y,curElem,lvl,p);
sigmah = sigma_h(x,y,curElem,lvl,p);

val = sum( (abs(sigmah - sigma0)).^(4/3),2 );
val = reshape(val,[1 1 length(x)]);

%% supply the error ||p_h-p||^2_L^(4/3) with p = DW(\grad u)
function val = L43errorPhminusP0Square(x,y,curElem,lvl,p)

%fprintf('SquareInit \n')
p.params.square=1;
sigma_0 = p.problem.sigma0;
sigma_h = p.statics.sigma_h;

sigma0 = sigma_0(x,y,curElem,lvl,p);
sigmah = sigma_h(x,y,curElem,lvl,p);

val = sum( (abs(sigmah - sigma0)).^(4/3),2 );
val = reshape(val,[1 1 length(x)]);

%% supply the error in the energy-functional \int_T |W(\grad u_h)+|u_h-f|^2-W(\grad u)-|u-f0|^2| 
function val = errorEnergy(x,y,curElem,lvl,p)

u_exact = p.problem.u_exact;
% I2u_h = p.statics.I2u_h;
gradU_exact = p.problem.gradU_exact;
% I2p_h = p.statics.I2p_h;
CONV = p.params.CONV;
if strcmp(CONV,'c')
nonLinearExact = p.problem.nonLinearReg;
else
nonLinearExact = p.problem.nonLinearExact;
end
f = p.problem.f0;

u_h = p.statics.u_h;
grad_h = p.statics.grad_h;

evalUexact = u_exact(x,y,p);
evalGradUexact = gradU_exact(x,y,p);

evalUh = u_h(x,y,curElem,lvl,p);
evalGradh = grad_h(x,y,curElem,lvl,p);

evalNonLinearExact = nonLinearExact(evalGradUexact(:,1),evalGradUexact(:,2),curElem,lvl,p);
evalNonLinearDiscrete = nonLinearExact(evalGradh(:,1),evalGradh(:,2),curElem,lvl,p);

% evalF = f(x,y,p);
evalF = f(x,y,curElem,lvl,p);

val = abs((evalNonLinearDiscrete - evalNonLinearExact) + ...
        abs(evalF-evalUh).^2-abs(evalF-evalUexact).^2);

val = reshape(val,[1 1 length(x)]);
