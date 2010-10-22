function p = drawDoF(dofNr,comp,localRes,h)

normRes = localRes*10;

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n4e = [1 2 3];
c4n = [0 0;
	   1 0;
	   0 1];
	   
c4n = h*c4n;
   
% c4n = [0 0;
% 		1 0;
% 		1/2 sqrt(2)/2];
p.level(1).geom.n4e = n4e;
p.level(1).geom.c4n = c4n;
p.level(1).geom.Db = [1 2];
p.level(1).geom.Nb = [];
p.PDE.lambda = 1;
p.PDE.mu = 1;

p = enumerate(p);
basisCoefficents = p.level(end).basisCoefficents;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;
hold on;

minX = 0; 
maxX = h;
minY = 0;
maxY = h;

XInorm = minX:h/normRes:maxX;
YInorm = minY:h/normRes:maxY;

XI = minX:h/localRes:maxX;
YI = minY:h/localRes:maxY;

curBasisCoefficents = basisCoefficents(:,:,1);


curSigma = zeros(24,1);
curSigma(dofNr) = 1;


curSigmaP3 = (curBasisCoefficents*curSigma)';


% UI = zeros(localRes);
% VI = zeros(localRes);
% for indX = 1:length(XI)
% 	for indY = 1:length(YI)
% 		val  = dofVal(XI(indX),YI(indY),comp,curSigmaP3,h);
% 		UI(indX,indY) = val(1);
% 		VI(indX,indY) = val(2);
% 	end
% end

UInorm = zeros(localRes);
VInorm = zeros(localRes);
for indX = 1:length(XInorm)
	for indY = 1:length(YInorm)
		val  = dofVal(XInorm(indX),YInorm(indY),comp,curSigmaP3,h);
		UInorm(indX,indY) = val(1)/sqrt(2);
		VInorm(indX,indY) = val(1)/sqrt(2);
	end
end

C = sqrt(UInorm.^2+VInorm.^2);
surf(XInorm,YInorm,zeros(size(C)),C,'EdgeColor','none');

% Draw Streamlines
%handle = streamslice(XI,YI,UI,VI);
%set(handle,'Color','black');

shading interp
% Plot Grid
%triplot(n4e,c4n(:,1),c4n(:,2),'LineWidth',1,'Color','k');

hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = dofVal(x,y,comp,coeffs,h)

%discrete sigma
if(x+y > h)
	val = [NaN,NaN];
else
	b1 = (h-x-y)/h;
	b2 = x/h;
	b3 = y/h;

	basis = [ b1;
			b2;
			b3;
			b1.*b2;
			b2.*b3;
			b1.*b3;
			b1.*b2.*(b1-b2);
			b2.*b3.*(b2-b3);
			-b1.*b3.*(b1-b3);	%%% TODO
			b1.*b2.*b3];


	I = [1:3:28];

	% coeff a^(j)
	coeffA = coeffs(I);
	% coeff b^(j)
	coeffB = coeffs(I+1);
	% coeff c^(j)
	coeffC = coeffs(I+2);


	sigmah = [coeffA*basis, coeffC*basis;...
			  coeffC*basis, coeffB*basis];

	%val = sigmah(comp,:);
% 	val = [norm(sigmah(1,:),2), norm(sigmah(2,:),2)];
	val = norm(sigmah);
end
