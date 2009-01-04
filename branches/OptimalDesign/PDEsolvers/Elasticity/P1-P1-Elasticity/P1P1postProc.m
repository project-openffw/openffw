function p = P1P1postProc(p)
% computes the post-processing datas for a P1-FE method 
% in linear elasticity. 

% Copyright 2007 Jan Reininghaus, David Guenther
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

%%%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%%
n4e = p.level(end).geom.n4e;

lambda = p.PDE.lambda;
mu = p.PDE.mu;

x = p.level(end).x;

nrNodes = p.level(end).nrNodes;
nrElems = p.level(end).nrElems;
grad4e = p.level(end).enum.grad4e;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = reshape(x(1:2*nrNodes)',nrNodes,2);

gradU = zeros(2,2,nrElems);
epsilon = zeros(2,2,nrElems);
sigma = zeros(2,2,nrElems);

genericUnit = eye(2);
u4e = zeros(3,2,nrElems);
for curElem = 1:nrElems
    curNodes = n4e(curElem,:);
    gradC = grad4e(:,:,curElem);
    curU = u(curNodes,:)';
	
	u4e(:,:,curElem) = curU';
    
    % calculation of piecewise constant gradient of u
    curGrad = curU*gradC;
    gradU(:,:,curElem) = curGrad;
	
    % calculation of piecewise constant strain (\eps(u))
    curStrain = (curGrad+curGrad')/2;
    epsilon(:,:,curElem) = curStrain;
    
    % calculation of piecewise constant stress (\sigma(u))
    traceUGrad = curGrad(1,1) + curGrad(2,2);
    curStress = lambda*traceUGrad*genericUnit + 2*mu*curStrain;
    sigma(:,:,curElem) = curStress;
end

%%%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%%%%%%%%%
p.level(end).u = u;
p.level(end).u4e = u4e;
p.level(end).gradU = gradU;
p.level(end).epsilon = epsilon;
p.level(end).sigma = sigma;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
