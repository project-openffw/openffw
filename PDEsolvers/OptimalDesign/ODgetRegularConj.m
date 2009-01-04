function p = ODgetRegularConj(p)

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
lambda = p.problem.lambda;
t1 = p.problem.t1;
t2 = p.problem.t2;

z = zeros(length(x),1);

amountMaterial = 0.5;

index1 = find( x <= (mu2+epsilon)*t1 );
z(index1) = 1/2*x(index1).^2/(mu2+epsilon) - amountMaterial*lambda*(mu1-mu2);

index2 = find( x > (mu2+epsilon)*t1 & x < (mu1+epsilon)*t2 );
z(index2) = mu2/2*t1^2 - amountMaterial*lambda*(mu1-mu2) + 1/2/epsilon*(mu2*t1-x(index2)).^2;

index3 = find( x >= (mu1+epsilon)*t2 );
z(index3) = 1/2*x(index3).^2/(mu1+epsilon) - mu1*t2/2*(t2-t1) - amountMaterial*lambda*(mu1-mu2);

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

index1 = find( x <= (mu2+epsilon)*t1 );
z(index1) = 1/(mu2+epsilon);

index2 = find( x > (mu2+epsilon)*t1 & x < (mu1+epsilon)*t2 );
z(index2) = 1/epsilon*(x(index2)-mu2*t1)./x(index2);

index3 = find( x >= (mu1+epsilon)*t2 );
z(index3) = 1/(mu1+epsilon);

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

index1 = find( x <= (mu2+epsilon)*t1 );
z(index1) = 1/(mu2+epsilon);

index2 = find( x > (mu2+epsilon)*t1 & x < (mu1+epsilon)*t2 );
z(index2) = 1/epsilon;

index3 = find( x >= (mu1+epsilon)*t2 );
z(index3) = 1/(mu1+epsilon);
