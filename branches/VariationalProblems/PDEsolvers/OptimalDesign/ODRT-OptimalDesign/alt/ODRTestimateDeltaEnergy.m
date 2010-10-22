function p = ODRTestimateDeltaEnergy(p)

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

h_T = max(length4ed(ed4e)')';
area4T = 0.25*h_T.^2; %works for triangles with 1 angle of 90° -> still to change;

% Lshape
% EnMid = -0.0963284556850688/3; %Aitken 2000, P1

% Lshape exact
% EnMid = -0.685115791830704/3; %Aitken 2000, P1
% EnMid = 1.80118931500645/3; %Aitken 700, f=1, P1

% SquareSlit
% EnMid = -0.146290348565174/4; %Aitken, P1

% SquareSlit exact
EnMid = 7.83238773919591/4; %Aitken 2000, P1

intEnergy = integrate(n4e,lvl,10,@getEnergy,p);
EnMid_h = sum(intEnergy);
fprintf('\nEnergy = %.15g\n',EnMid_h)
nu = intEnergy - EnMid.*area4T;

%% OUTPUT
p.level(end).etaT = nu;
p.level(end).estimatedError = abs(sum(nu));%norm(nu,2);

%% supply the discrete energy
function val = getEnergy(x,y,curElem,lvl,p)

energy_h = p.statics.energy_h;
val = energy_h(x,y,curElem,lvl,p);
