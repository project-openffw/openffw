function p = getNonLinearReg(p)

p.problem.nonLinearReg = @nonLinearReg;
p.problem.nonLinearRegDerA = @nonLinearRegDerA;
p.problem.nonLinearRegDerB = @nonLinearRegDerB;
p.problem.nonLinearRegSecDerA = @nonLinearRegSecDerA;
p.problem.nonLinearRegSecDerB = @nonLinearRegSecDerB;
p.problem.nonLinearRegSecDerC = @nonLinearRegSecDerC;

%% exact nonlinear operator (relaxed)
%% W**(F) = (max{|F|^2-1,0})^2 + 4(|F|^2 - (F2,F)^2)
function z = nonLinearReg(x,y,curElem,lvl,p)

F1 = p.problem.F1;
F2 = p.problem.F2;
z = zeros(length(x),1);
maxvalue = max(0,x.^2+y.^2-1);
z = maxvalue.^2 + 4*(x.^2+y.^2 - (F2(1)*x + F2(2)*y).^2);

%% exact derivative of nonlinear operator (relaxed)
%% DW**(F)[G] = W'**1(F)(F,G) + W'**2(F)(F2,F)(F2,G)
%% W'**1(F) = 4 max{|F|^2-1,0} + 8
function z = nonLinearRegDerA(x,y,curElem,lvl,p)

z = zeros(length(x),1);
maxvalue = max(0,x.^2+y.^2-1);
z = 4*maxvalue + 8*ones(length(x),1);

%% W'**2(F) = -8
function z = nonLinearRegDerB(x,y,curElem,lvl,p)

z = -8*ones(length(x),1);

%% exact second derivative of nonlinear operator (relaxed)
%% D2W**(F)[G,H] = W''**1(F)(G,H) + W''**2(F)(F,G)(F,H) +
%%                 W''**3(F)(F2,H)(F2,G)
%% W''**1(F) = 4 max{|F|^2-1,0} + 8
function z = nonLinearRegSecDerA(x,y,curElem,lvl,p)

maxvalue = max(0,x.^2+y.^2-1);
z = 4*maxvalue + 8*ones(length(x),1);

%% W''**2(F) = 8*HF(|F|^2-1)
%% HF(x)=0 for x<=0 and HF(x)=1 for x > 0, Heavyside-function
function z = nonLinearRegSecDerB(x,y,curElem,lvl,p)

z = zeros(length(x),1);

val = x.^2+y.^2-1;
Hi1 = find(val <= 0);
z(Hi1) = 0;
Hi2 = find(val>0);
z(Hi2) = 8;

%% W''**3(F) = -8
function z = nonLinearRegSecDerC(x,y,curElem,lvl,p)

z = -8*ones(length(x),1);












