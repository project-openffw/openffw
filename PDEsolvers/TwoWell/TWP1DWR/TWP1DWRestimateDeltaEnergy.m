function p = TWP1estimateDeltaEnergy(p)

% Copyright 2008 Lena Noack
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
% length4ed = p.level(end).enum.length4ed;
lvl = size(p.level,2);
n4e = p.level(lvl).geom.n4e;
% n4ed = p.level(end).enum.n4ed;
% ed4e = p.level(end).enum.ed4e;
curLvl = length(p.level);
degree = 1;
nrNodes = p.level(end).nrNodes;
length4ed = p.level(end).enum.length4ed;
ed4e = p.level(end).enum.ed4e;
problem = p.params.problem.name;
CONV = p.params.CONV;

h_T = max(length4ed(ed4e)')';
area4T = 0.25*h_T.^2; %works for triangles with 1 angle of 90° -> still to change;

if strcmp(problem,'TwoWell_Square') % TWCF
  if strcmp(CONV,'c')
    EnMid = 1.45193218965457/1.5; %Aitken 1400, Convex, f0=1, uD=0
  else
    EnMid = 1.67112143088421; %Aitken 1400, Nonconvex, f0=1, uD=0
  end
else % TWCO
  if strcmp(CONV,'c')
    EnMid = 0.107826374840891/1.5; %Aitken 1400, Convex, f0, uD given
  else
    EnMid = 1.07002831984151/1.5; %Aitken 1400, Nonconvex, f0, uD given
  end
end

intEnergy = integrate(n4e,lvl,10,@getEnergy,p);
EnMid_h = sum(intEnergy);
nu = intEnergy - EnMid.*area4T;
fprintf('\nEnergy = %.15g\n',EnMid_h)

%% OUTPUT
p.level(end).etaT = nu;
%p.level(end).estimatedError = abs(sum(nu));%norm(nu,2);        %est. for ||sigma-sigma_h||^2
p.level(end).estimatedError = (abs(sum(nu))).^0.5;%norm(nu,2);  %est. for ||sigma-sigma_h||

%% supply the discrete energy
function val = getEnergy(x,y,curElem,lvl,p)

energy_h = p.statics.energy_h;
val = energy_h(x,y,curElem,lvl,p);
