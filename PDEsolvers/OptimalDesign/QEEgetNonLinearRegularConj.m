function p = getNonLinearRegularConj(p)

p.problem.nonLinearFunc = @nonLinearFunc;
p.problem.nonLinearFuncDer = @nonLinearFuncDer;
p.problem.nonLinearFuncSecDer = @nonLinearFuncSecDer;

%% W_\epsilon(x)
function z = nonLinearFunc(x,curElem,lvl,p)

z = zeros(length(x),1);
z = x.*x;

%% DW_\epsilon(x)/x
function z = nonLinearFuncDer(x,curElem,lvl,p)

z = ones(length(x),1);

%% D^2W_\epsilon
function z = nonLinearFuncSecDer(x,curElem,lvl,p)

z = ones(length(x),1);
