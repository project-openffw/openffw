function [x,residual,jacobian,t] = ODlineSearch(x,lvl,p)
% author: David Guenther 

pdeSolver = p.params.pdeSolver;
getNonLinearSolution = str2func([pdeSolver,'getNonLinearSolution']);

[residual, jacobian] = getNonLinearSolution(x,p);

tolerance = loadField('p.params','tolerance',p,1e-6);

eta = norm(residual);

fprintf('\nNewton-Raphson Scheme: \n');
fprintf('====================== \n');
fprintf('%i, %1.12f \n',0,eta);

t = 1;
nrInterations = 5;

for k = 1:nrInterations
    if k > 1
        [residual, jacobian] = getNonLinearSolution(x,p);
    end
    correction = jacobian \ residual;
    [x, residual, t] = lineSearchStepSize(x,correction,eta(end),t,p);
    eta = [eta, norm(residual)];
    fprintf('%i, %1.12f \n',size(eta,2)-1,eta(end));
    if eta(end) < tolerance;
        return
    end
end




function [xNew,residual,t] = lineSearchStepSize(x,correction,etaOld,t,p)

pdeSolver = p.params.pdeSolver;
getNonLinearSolution = str2func([pdeSolver,'getNonLinearSolution']);

t = 2*t;
eta = 1e9;

while eta > etaOld
    t = t/2;
    xNew = x - t*correction;
    [residual, dontUse] = getNonLinearSolution(xNew,p);
    eta = norm(residual);
end

