function p = getConj2(p)

p.problem.conjNonLinearFunc = @conjNonLinearFunc;
p.problem.conjNonLinearFuncDer = @conjNonLinearFuncDer;
p.problem.conjNonLinearFuncSecDer = @conjNonLinearFuncSecDer;

%% W***(x)
function z = conjNonLinearFunc(x,y,curElem,lvl,p)
z = ones(length(x),1);

%% exact derivative of conjugative nonlinear operator (relaxed)
%% DW***(F)[G] = (W***'(F),G)
function [z1 z2] = conjNonLinearFuncDer(x,y,curElem,lvl,p)
z1 = ones(length(x),1);
z2 = ones(length(x),1);

%% exact derivative of conjugative nonlinear operator (relaxed)
%% D2W***(F)[G,H] = W***''(F) * (G,H)
function z = conjNonLinearFuncSecDer(x,y,curElem,lvl,p)
z = ones(length(x),1);
