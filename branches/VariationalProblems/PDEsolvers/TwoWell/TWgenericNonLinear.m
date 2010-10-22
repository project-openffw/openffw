function p = genericNonLinear(p)

p.problem.nonLinearExact = @nonLinearExact;
p.problem.nonLinearExactDerA = @nonLinearExactDerA;
p.problem.nonLinearExactDerB = @nonLinearExactDerB;
p.problem.nonLinearExactSecDerA = @nonLinearExactSecDerA;
p.problem.nonLinearExactSecDerB = @nonLinearExactSecDerB;

%% exact nonlinear operator (not relaxed)
%% W(F) = |F-F1|^2 |F-F2|^2
function z = nonLinearExact(x,y,curElem,lvl,p)

F1 = p.problem.F1;
F2 = p.problem.F2;
z = zeros(length(x),1);
z = ((x-F1(1)).^2 + (y-F1(2)).^2).*((x-F2(1)).^2 + (y-F2(2)).^2);

%% exact derivative of nonlinear operator (not relaxed)
%% DW(F)[G] = W'1(F)(F-F1,G) + W'2(F)(F-F2,G)
%% W'1(F) = 2|F-F2|^2
function z = nonLinearExactDerA(x,y,curElem,lvl,p)

F1 = p.problem.F1;
F2 = p.problem.F2;
z = zeros(length(x),1);
z = 2*((x-F2(1)).^2 + (y-F2(2)).^2);

%% W'2(F) = 2|F-F1|^2
function z = nonLinearExactDerB(x,y,curElem,lvl,p)

F1 = p.problem.F1;
F2 = p.problem.F2;
z = zeros(length(x),1);
z = 2*((x-F1(1)).^2 + (y-F1(2)).^2);

%% exact second derivative of nonlinear operator (not relaxed)
%% D2W(F)[G,H] = W''1(F)(G,H) + W''2(F)((F-F1,G)*(F-F2,H) + (F-F1,H)*(F-F2,G))
%% W''1(F) = 2|F-F1|^2 + 2|F-F2|^2
function z = nonLinearExactSecDerA(x,y,curElem,lvl,p)

F1 = p.problem.F1;
F2 = p.problem.F2;
z = zeros(length(x),1);
z = 2*((x-F1(1)).^2 + (y-F1(2)).^2) + 2*((x-F2(1)).^2 + (y-F2(2)).^2);

%% W''2(F) = 4
function z = nonLinearExactSecDerB(x,y,curElem,lvl,p)

z = 4*ones(length(x),1);













