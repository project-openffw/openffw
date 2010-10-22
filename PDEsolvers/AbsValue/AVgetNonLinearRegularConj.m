function p = getNonLinearRegularConj(p)

%p.problem.nonLinearFunc = @nonLinearFunc;
%p.problem.nonLinearFuncDer = @nonLinearFuncDer;
%p.problem.nonLinearFuncSecDer = @nonLinearFuncSecDer;

%% W_\epsilon(x)
function z = nonLinearFunc(x,curElem,lvl,p)

solver = p.params.solver;

if strcmp(solver,'eps-h-dependence')
    area4e = p.level(lvl).enum.area4e;
    epsPower = p.params.epsPower;
    epsilon = area4e(curElem)^epsPower;
else
    epsilon = p.problem.epsilon;
end

z = zeros(length(x),1);

index1 = find(x < 2);
z(index1) = (1+2*epsilon)*x(index1).^2/4;
index2 = find(x >= 2 & x <= 4);
z(index2) = epsilon/2*x(index2).^2 + x(index2) - 1;
index3 = find(x > 4);
z(index3) = (1+4*epsilon)*x(index3).^2/8 + 1;

%% DW_\epsilon(x)/x
function z = nonLinearFuncDer(x,curElem,lvl,p)

solver = p.params.solver;

if strcmp(solver,'eps-h-dependence')
    area4e = p.level(lvl).enum.area4e;
    epsPower = p.params.epsPower;
    epsilon = area4e(curElem)^epsPower;
else
    epsilon = p.problem.epsilon;
end

z = zeros(length(x),1);

index1 = find(x < 2);
z(index1) = (1+2*epsilon)/2;
index2 = find(x >= 2 & x <= 4);
z(index2) = (epsilon*x(index2) + 1)./x(index2);
index3 = find(x > 4);
z(index3) = (1+4*epsilon)/4;

%% D^2W_\epsilon
function z = nonLinearFuncSecDer(x,curElem,lvl,p)

solver = p.params.solver;

if strcmp(solver,'eps-h-dependence')
    area4e = p.level(lvl).enum.area4e;
    epsPower = p.params.epsPower;
    epsilon = area4e(curElem)^epsPower;
else
    epsilon = p.problem.epsilon;
end

z = zeros(length(x),1);

index1 = find(x < 2);
z(index1) = (1+2*epsilon)/2;
index2 = find(x >= 2 & x <= 4);
z(index2) = epsilon;
index3 = find(x > 4);
z(index3) = (1+4*epsilon)/4;
