function p = genericNonLinear(p)

%p.problem.nonLinearConjExact = @nonLinearConjExact;
p.problem.nonLinearExact = @nonLinearExact;
p.problem.nonLinearExactDer = @nonLinearExactDer;
p.problem.nonLinearExactSecDer = @nonLinearExactSecDer;

%% exact conjugate nonlinear operator W*(x) (not regularised)
function z = nonLinearConjExact(x,curElem,lvl,p)

z = 0;

%% exact nonlinear operator W(x) (not regularised)
%% W(|Du|)=|Du|^p + regParam*|Du|^2
function z = nonLinearExact(x,curElem,lvl,p)

regParam = p.problem.regParam;
pv = p.problem.pv;
z = zeros(length(x),1);
z = x.^pv + regParam*x.^2;

%% exact derivative of nonlinear operator, W'(x)/x (not regularised)
%% DW(|Du|)[Dv] = W'(|Du|) * (Du,Dv)_R^n
%% W'(|Du|) = p*|Du|^(p-2) + regParam*2
function z = nonLinearExactDer(x,curElem,lvl,p)

regParam = p.problem.regParam;
pv = p.problem.pv;
z = pv*x.^(pv-2) + regParam*2;

%% exact second derivative of nonlinear operator, W''(x) (not regularised)
%% D^2W(|Du|)[Dv,Dw] = W''(|Du|)*(Dv,Du)_R^n*(Dw,Du)_R^n +
%% (W'(|Du|)/|Du|)*(Dv,Dw)_R^n
%% W''(|Du|) = p(p-2)*|Du|^(p-4) + regParam*2
function z = nonLinearExactSecDer(x,curElem,lvl,p)

regParam = p.problem.regParam;
pv = p.problem.pv;
if pv == 2
    z = zeros(length(x),1);
else
    z = pv*(pv-2)*x.^(pv-4) + regParam*2;
end
