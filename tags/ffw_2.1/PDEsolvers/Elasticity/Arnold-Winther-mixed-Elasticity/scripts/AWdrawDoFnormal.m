function p = AWdrawDoFnormal(dofNr,normal,localRes)
% draw DoF for normal

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


normRes = localRes*10;

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n4e = [1 2 3];
% c4n = [0 0;
% 	   1 0;
% 	   0 1];
c4n = [-1/2 0;
	   1/2 0;
	   0	sqrt(3)/2];
p.level(1).geom.n4e = n4e;
p.level(1).geom.c4n = c4n;
p.level(1).geom.Db = [1 2];
p.level(1).geom.Nb = [];
p.PDE.lambda = 1;
p.PDE.mu = 1;

p = AWenumerate(p);
basisCoefficents = p.level(end).basisCoefficents;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;
hold on;

minX = -1/2; 
maxX = 1/2;
minY = 0;
maxY = sqrt(3)/2;

XInorm = minX:1/normRes:maxX;
YInorm = minY:1/normRes:maxY;


curBasisCoefficents = basisCoefficents(:,:,1);


curSigma = zeros(24,1);
curSigma(dofNr) = 1;


curSigmaP3 = (curBasisCoefficents*curSigma)';



UVnorm = zeros(localRes);
for indX = 1:length(XInorm)
	for indY = 1:length(YInorm)
		val  = dofVal(XInorm(indX),YInorm(indY),normal,curSigmaP3);
		UVnorm(indX,indY) = norm(val,2);
	end
end

C = UVnorm;
surf(XInorm,YInorm,zeros(size(C')),C','EdgeColor','none');

% Draw Streamlines
% handle = streamslice(XI,YI,UI,VI);
% set(handle,'Color','black');

shading interp
% Plot Grid
triplot(n4e,c4n(:,1),c4n(:,2),'LineWidth',1,'Color','k');

hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = dofVal(x,y,normal,coeffs)
coords = [-1/2 0;
	   1/2 0;
	   0	sqrt(3)/2];

% normals = [0 -1;
% 		   sqrt(2)/2, sqrt(2)/2;
% 		   -1 0]';
	   
normals = [0 -1;
		   sqrt(2)/2, sqrt(2)/2;
		   -1 0]';
	   
curNormal = normals(:,normal);
%discrete sigma
% if(x+y > 1)
if(sqrt(3)*x+y > sqrt(3)/2 || y-sqrt(3)*x > sqrt(3)/2)
	val = [NaN,NaN];
else
% 	P1 = coords(1,:);
% M = [coords(2,:) - coords(1,:)
% 	 coords(3,:) - coords(1,:)];
% 	z = inv(M)*([x y] - P1)';
% 	x = z(1);
% 	y = z(2);
% 	b1 = 1-x-y;
% 	b2 = x;
% 	b3 = y;
	
	b1 = -x - y/sqrt(3) + 1/2;
	b2 = x - y/sqrt(3) + 1/2;
	b3 = 2*y/sqrt(3);

	basis = [ b1;
			b2;
			b3;
			b1.*b2;
			b2.*b3;
			b1.*b3;
			b1.*b2.*(b1-b2);
			b2.*b3.*(b2-b3);
			b1.*b3.*(b3-b1);
			b1.*b2.*b3];


% 	I = [1:3:28];
% 
% 	% coeff a^(j)
% 	coeffA = coeffs(I);
% 	% coeff b^(j)
% 	coeffB = coeffs(I+1);
% 	% coeff c^(j)
% 	coeffC = coeffs(I+2);

	I = [1:10];

	% coeff a^(j)
	coeffA = coeffs(I);
	% coeff b^(j)
	coeffB = coeffs(I+10);
	% coeff c^(j)
	coeffC = coeffs(I+20);


	sigmah = [coeffA*basis, coeffC*basis;...
			  coeffC*basis, coeffB*basis];

% 	val = sigmah*curNormal;
	val = norm(sigmah,2);
end
