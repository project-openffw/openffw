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

% SquareSlit
EnMid = -0.147276354410246/4; %Aitken 3200 P0-RT0
sigmaSqL2A = 0.330929232065413; %Aitken 3200 P0-RT0

% SquareSlit exact
%EnMid = 0.1/4; %Aitken 3200 P0-RT0
%sigmaSqL2A = 0.1; %Aitken 3200 P0-RT0

% Lshape
%EnMid = -0.0963349567479233/3; %Aitken 3200 P0-RT0
%sigmaSqL2A = 0.216858427951399; %Aitken 3200 P0-RT0

% Lshape exact
%EnMid = 0.1/3; %Aitken 3200 P0-RT0
%EnMid = 0.1/3; %Aitken 3200 P0-RT0, f=1
%sigmaSqL2A = 0.1; %Aitken 3200 P0-RT0

intEnergy = integrate(n4e,lvl,10,@getEnergy,p);
EnMid_h = sum(intEnergy);
fprintf('\nEnergy = %.15g\n',EnMid_h)
nu = intEnergy - EnMid.*area4T;

sigmaSquare = integrate(n4e,lvl,10,@getSigma,p);
sigmaSqL2 = sum(sigmaSquare);                      %||sigma_h||^2_L^2(Omega)
fprintf('\nSigma = %.15g\n',sigmaSqL2)

%% OUTPUT
p.level(end).sigmaSqL2 = sigmaSqL2;
p.params.sigmaSqL2A = sigmaSqL2A;
p.level(end).etaT = nu;
%p.level(end).estimatedError = abs(sum(nu));%norm(nu,2);       %est. for ||sigma-sigma_h||^2
p.level(end).estimatedError = (abs(sum(nu))).^0.5;%norm(nu,2); %est. for ||sigma-sigma_h||

%% supply the discrete energy
function val = getEnergy(x,y,curElem,lvl,p)

energy_h = p.statics.energy_h;
val = energy_h(x,y,curElem,lvl,p);

%% supply the discrete energy
function val = getSigma(x,y,curElem,lvl,p)

sigma_h = p.statics.sigma_h;
sigmah = sigma_h(x,y,curElem,lvl,p);

sigmah = reshape(sigmah',[2 1 length(x)]);

val = sum(sigmah.^2,1);