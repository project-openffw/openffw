function p = getRegularConj(p)

    p.problem.conjNonLinearFunc = @conjNonLinearFunc;
    p.problem.conjNonLinearFuncDer = @conjNonLinearFuncDer;
    p.problem.conjNonLinearFuncSecDer = @conjNonLinearFuncSecDer;

%% W*
function z = conjNonLinearFunc(x,curElem,lvl,p)

pv = p.problem.pv;
if pv==2
    z = -0.5*x.^2;
else
    fprintf('Conjugate Functions can only be calculated for p=2 (and \varepsilon=0)!');
end

%% DW*
function z = conjNonLinearFuncDer(x,curElem,lvl,p)

z = -x;

%% D^2W*
function z = conjNonLinearFuncSecDer(x,curElem,lvl,p)

z = -ones(length(x),1);
