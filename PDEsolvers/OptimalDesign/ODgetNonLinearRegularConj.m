function p = getNonLinearRegularConj(p)
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

p.problem.nonLinearFunc = @nonLinearFunc;
p.problem.nonLinearFuncDer = @nonLinearFuncDer;
p.problem.nonLinearFuncSecDer = @nonLinearFuncSecDer;

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
