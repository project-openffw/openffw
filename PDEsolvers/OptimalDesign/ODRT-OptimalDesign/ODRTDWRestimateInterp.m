function p = ODRTestimateInterp(p)

% Copyright 2008 Joscha Gedicke, Lena Noack
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

%% INPUT
length4ed = p.level(end).enum.length4ed;
n4e = p.level(end).geom.n4e;
n4ed = p.level(end).enum.n4ed;
ed4e = p.level(end).enum.ed4e;
curLvl = length(p.level);
degree = 1;
nrNodes = p.level(end).nrNodes;


%% ESTIMATE
h_T = max(length4ed(ed4e),[],2);

nu_T = integrate(n4e,curLvl,2*degree,@Residuum,p);
nu_E = integrate(n4ed,curLvl,degree,@NormalJumpSigma,p);
primalResiduum = sqrt( nu_T + 1/4*1./h_T.*sum(nu_E(ed4e),2) );

p = P2RTInterpolationDual(p);
dnu_T = integrate(n4e,curLvl,4*degree,@funcHandleDualWeightElems,p);
dnu_E = integrate(n4ed,curLvl,4*degree,@funcHandleDualWeightEdges,p);
dualWeight = sqrt( dnu_T + h_T.*sum(dnu_E(ed4e),2) );


nu = primalResiduum.*dualWeight;

%estErrorT = integrate(n4e,curLvl,4*degree,@funcHandleWeightsEstError,p)
%estErrorE = integrate(n4ed,curLvl,4*degree,@funcHandleWeightsEstErrorEd,p)

%% OUTPUT
p.level(end).etaT = nu;
p.level(end).estimatedError = sqrt(sum(nu));  %est. for ||sigma-sigma_h||
%p.level(end).estimatedError = abs(sum(estimatedError));


function val = Residuum(x,y,curElem,curLvl,p)
f = p.problem.f;
curf     = f(x,y,curElem,curLvl,p);

residuum = - curf(:);
val(1,:,:) = (residuum.^2)';

function val = NormalJumpSigma(x,y,curEdge,lvl,p)
sigma_h = p.statics.sigma_h;
e4ed = p.level(lvl).enum.e4ed;
normals4ed = p.level(lvl).enum.normals4ed;

elems = e4ed(curEdge,:);
normal = normals4ed(curEdge,:);

evalSigma1 = sigma_h(x,y,elems(1),lvl,p)*normal';

if elems(2) ~= 0
    % E is an interior edge
    evalSigma2 = sigma_h(x,y,elems(2),lvl,p)*normal';
else
    % E is a Dirichlet Edge
    evalSigma2 = evalSigma1;
end

val = zeros(1,1,length(x));
val(1,1,:) = sum((evalSigma1 - evalSigma2).^2,2);


