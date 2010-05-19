function p = ODRTestimate_Jump(p)
%% INPUT
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
n4ed = p.level(end).enum.n4ed;
length4ed = p.level(end).enum.length4ed;

lvl = size(p.level,2);

degree = loadField('p.params','nonLinearExactIntegrateDegree',p,19);
%% estimation of the error
eta4ed = length4ed.*integrate(n4ed,lvl,degree,@intTangentJump,p);

%% compose the error terms
estimatedError = sqrt( sum(eta4ed) );

%% OUTPUT
p.level(end).etaEd = sqrt(eta4ed);
p.level(end).estimatedError = estimatedError;

%% supply tangent jump: int_E ([DW*(p)]*tau)^2
function val = intTangentJump(points,curEdge,lvl,p)

grad_h = p.statics.grad_h;
e4ed = p.level(lvl).enum.e4ed;
tangents4ed = p.level(lvl).enum.tangents4ed;

elems = e4ed(curEdge,:);
tangent = tangents4ed(curEdge,:)';

evalGrad_1 = grad_h(points,elems(1),lvl,p);
if elems(2) ~= 0
    evalGrad_2 = grad_h(points,elems(2),lvl,p);
else
    evalGrad_2 = evalGrad_1;
end

diffGrad = (evalGrad_1 - evalGrad_2)*tangent;
val = diffGrad.^2;
val = reshape(val,[1 1 length(points(:,1))]);
