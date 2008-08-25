function p = CRestimate(p)
% estimate the energy error for a nonconforming CR-FE method 
% through computing the oscillations of f, and
% the jump of the flux. 

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


%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load discrete solution
ph = p.level(end).grad4e;

% problem data
f = p.problem.f;
g = p.problem.g;

% load enumerated data
e4ed = p.level(end).enum.e4ed;
NbEd = p.level(end).enum.NbEd;
length4ed = p.level(end).enum.length4ed;
area4e = p.level(end).enum.area4e;
midpoint4e = p.level(end).enum.midPoint4e;
midPoint4ed = p.level(end).enum.midPoint4ed;
normals4NbEd = p.level(end).enum.normals4NbEd;

% load integration parameters
degreeJumpTerm = p.params.integrationDegrees.estimate.jumpTerm;
degreeVolumeTerm = p.params.integrationDegrees.estimate.volumeTerm;
degreeOscTerm = p.params.integrationDegrees.estimate.oscTerm;

%% Estimate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% jump term
phU = ph(:,1);
phV = ph(:,2);

e4ed = sort(e4ed,2);
[I,dontUse] = find(e4ed == 0);
e4ed(I,1) = e4ed(I,2);

ph4edU1 = phU(e4ed(:,1));
ph4edU2 = phU(e4ed(:,2));
ph4edU = ph4edU1 - ph4edU2;

ph4edV1 = phV(e4ed(:,1));
ph4edV2 = phV(e4ed(:,2));
ph4edV = ph4edV1 - ph4edV2;

jump4ed = length4ed.^2 .* (ph4edU.^2 + ph4edV.^2);

% Neumann boundary
if ~isempty(NbEd) 
	g4NbEd = g(midPoint4ed(NbEd,:),normals4NbEd,p);
	ph4Nb = [ph4edU1(NbEd),ph4edV1(NbEd)];
	ph4NbNormal = sum(ph4Nb .* normals4NbEd,2);
	
	jump4ed(NbEd) = length4ed(NbEd).^2 .* (g4NbEd - ph4NbNormal).^2;
end


% Oscillations 
% boundary elements are considered to be their own neighbours
[I,dontUse] = find(e4ed == 0);
e4ed(I,2) = e4ed(I,1);

firstElem = e4ed(:,1);
secElem = e4ed(:,2);

area4firstElem = area4e(firstElem);
area4secElem = area4e(secElem);

f4firstElem = f(midpoint4e(firstElem,:),p);
f4secElem = f(midpoint4e(secElem,:),p);

f4patch = (area4firstElem.*f4firstElem + area4secElem.*f4secElem)./...
            (area4secElem + area4secElem);

f4T1 = area4e(firstElem).*(f4firstElem - f4patch).^2;
f4T2 = area4e(secElem).*(f4secElem - f4patch).^2;

osc4ed = length4ed.^2 .* ( f4T1 + f4T2 );


% estimate
eta4ed = sqrt(jump4ed + osc4ed);

%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).jump4ed = sqrt(jump4ed);
p.level(end).osc4ed = sqrt(osc4ed);
p.level(end).etaEd = eta4ed;
p.level(end).estimatedError = norm(eta4ed,2);