function p = QEEgetRegularConj(p)

p.problem.conjNonLinearFunc = @conjNonLinearFunc;
p.problem.conjNonLinearFuncDer = @conjNonLinearFuncDer;
p.problem.conjNonLinearFuncSecDer = @conjNonLinearFuncSecDer;

%% W*_\epsilon(x)
function z = conjNonLinearFunc(x,curElem,lvl,p)

solver = p.params.solver;

if strcmp(solver,'ODeps_h_dependence')
    area4e = p.level(lvl).enum.area4e;
    epsPower = p.params.epsPower;
    epsFactor = p.params.epsFactor;
    epsilon = epsFactor*area4e(curElem)^(epsPower/2);
else
    epsilon = p.problem.epsilon;
end

mu1 = p.problem.mu1;
mu2 = p.problem.mu2;
t1 = p.problem.t1;
t2 = p.problem.t2;

z = zeros(length(x),1);

k = 1/epsilon * (1/mu1 + 1/epsilon)^(-1);

index1 = find( x <= mu2/k *t1 );
z(index1) = 1/(2*mu2)* x(index1).^2*k^2 + 1/(2*epsilon) * x(index1) * abs(1-k);

index2 = find( x > mu2/k *t1 );
z(index2) = 1/(2*mu1)* x(index2).^2*k^2 + 1/(2*epsilon) * x(index2) * abs(1-k)...
            +1/2 * mu2 * t1^2 - 1/2 * mu1 * t2^2;

%% DW*_\epsilon(x)/x
function z = conjNonLinearFuncDer(x,curElem,lvl,p)

solver = p.params.solver;

if strcmp(solver,'ODeps_h_dependence')
    area4e = p.level(lvl).enum.area4e;
    epsPower = p.params.epsPower;
    epsFactor = p.params.epsFactor;
    epsilon = epsFactor*area4e(curElem)^(epsPower/2);
else
    epsilon = p.problem.epsilon;
end

mu1 = p.problem.mu1;
mu2 = p.problem.mu2;
t1 = p.problem.t1;
t2 = p.problem.t2;

z = zeros(length(x),1);

k = 1/epsilon * (1/mu1 + 1/epsilon)^(-1);

index1 = find( x <= mu2/k *t1 );
z(index1) = k^2/(mu2)+ 1/(2*epsilon* x(index1)) * abs(1-k);

index2 = find( x > mu2/k *t1 );
z(index2) = k^2/(mu1)+ 1/(2*epsilon* x(index2)) * abs(1-k);

%% D^2W*_\epsilon
function z = conjNonLinearFuncSecDer(x,curElem,lvl,p)

solver = p.params.solver;

if strcmp(solver,'ODeps_h_dependence')
    area4e = p.level(lvl).enum.area4e;
    epsPower = p.params.epsPower;
    epsFactor = p.params.epsFactor;
    epsilon = epsFactor*area4e(curElem)^(epsPower/2);
else
    epsilon = p.problem.epsilon;
end

mu1 = p.problem.mu1;
mu2 = p.problem.mu2;
t1 = p.problem.t1;
t2 = p.problem.t2;

z = zeros(length(x),1);

k = 1/epsilon * (1/mu1 + 1/epsilon)^(-1);

index1 = find( x <= mu2/k *t1 );
z(index1) = k^2/(mu2);

index2 = find( x > mu2/k *t1 );
z(index2) = k^2/(mu1);
