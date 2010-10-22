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

h_T = max(length4ed(ed4e)')';
area4T = 0.25*h_T.^2; %works for triangles with 1 angle of 90° -> still to change;

if strcmp(problem,'AbsValue_Square') % Square
    EnMid = -0.00878618555043987; %Aitken 900, p=2, regParam=0.0, |Omega|=1 
    %EnMid = -0.00798744140949077; %Aitken 900, p=2, regParam=0.1, |Omega|=1 
elseif strcmp(problem,'AbsValue_Square_exact') % Square exact
    EnMid = -0.0244777167125984; %Aitken 3900, p=2, regParam=0.1, |Omega|=1
elseif strcmp(problem,'AbsValue_Lshape') % Lshape
    EnMid = -0.0535018557045206/3; %Aitken 2900, p=2, regParam=0.0, |Omega|=3
    %EnMid = -0.0486292609418941/3; %Aitken 700, p=2, regParam=0.1, |Omega|=3
else % Lshape exact
    EnMid = 1.94258379192296/3; %Aitken 2900, p=2, regParam=0.1, |Omega|=3
end

intEnergy = integrate(n4e,lvl,10,@getEnergy,p);
EnMid_h = sum(intEnergy);
fprintf('\nEnergy = %.15g\n',EnMid_h)
nu = intEnergy - EnMid.*area4T;

%% OUTPUT
p.level(end).etaT = nu;
%p.level(end).estimatedError = abs(sum(nu));%norm(nu,2); %est. for ||sigma-sigma_h||^2
p.level(end).estimatedError = (abs(sum(nu))).^(1/2);     %est. for ||sigma-sigma_h||

%% supply the discrete energy
function val = getEnergy(x,y,curElem,lvl,p)

energy_h = p.statics.energy_h;
val = energy_h(x,y,curElem,lvl,p);
