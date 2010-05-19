function [x,residual,jacobian,t] = ODlineSearch(x,lvl,p)
% author: David Guenther 
% Copyright 2007 David Guenther
%
% This file is part of FFW.
%
% FFW is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% FFW is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
%

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

