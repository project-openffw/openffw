function evalExamples2(p,k,refineType,linestyle)

if nargin < 3
    error('the input arguments are: the datastructure p, the number (integer) of the geometrie (1: square, 2: l-shape, 3: octagon, 4: slitted square), and the refinement type (string) (uniform, adaptive).')
end

if nargin < 4
    linestyle = '-';
end

if k == 1
    % square
    name = 'Square';
    limit = 0.027983242756785;
elseif k == 2
    % lshape
    name = 'Lshape';
    limit = 0.160674590782529;
elseif k == 3
    % octagon
    name = 'Octagon';
    limit = 0.257317403528310;
elseif k == 4
    % square slit
    name = 'Square Slit';
    limit = 0.247295558950406;
end

nrLevels = length(p.level);

areaContactZone = zeros(nrLevels-1,1);
errorAverage = zeros(nrLevels-1,1);
errorProj = zeros(nrLevels-1,1);
errorJump = zeros(nrLevels-1,1);
errorEnergy = zeros(nrLevels-1,1);
nrDoF = zeros(nrLevels-1,1);

for lvl = 2:nrLevels
    areaContactZone(lvl-1) = p.level(lvl).areaContactZone;
    errorAverage(lvl-1) = p.level(lvl).phminusAph;
    errorProj(lvl-1) = p.level(lvl).estimatedError;
%     errorJump(lvl-1) = p.level(lvl).errorJump;
    errorEnergy(lvl-1) = p.level(lvl).conjEnergy - limit;
    nrDoF(lvl-1) = p.level(lvl).nrDoF;
end

firstIndex = find(nrDoF > 500,1,'first');

index = firstIndex:nrLevels-1;
nrDoF = nrDoF(index);

term1 = errorEnergy(index).^(1/2);
term2 = areaContactZone(index).^(1/2);
term3 = errorAverage(index);
term4 = ((areaContactZone(index) + errorAverage(index)).*errorAverage(index)).^(1/2);
term5 = errorProj(index);
term6 = errorJump(index);
matrix = [index',nrDoF,term1,term2,term3,term4,term5,term6];
latex(matrix)

hold all
% loglog(nrDoF,term1,'Displayname',[refineType,' ','(E(\sigma_{\theta,h})-E(\sigma))^{1/2}'],'Marker','s','Markersize',10,'Linewidth',2,'Linestyle',linestyle)
loglog(nrDoF,term2,'Displayname',[refineType,' ','Area Contact Zone'],'Marker','d','Markersize',10,'Linewidth',2,'Linestyle',linestyle)
loglog(nrDoF,term3,'Displayname',[refineType,' ','||\sigma_{h,\theta}-A\sigma_{h,\theta}||_{L^2}'],'Marker','p','Markersize',10,'Linewidth',2,'Linestyle',linestyle)
loglog(nrDoF,term4,'Displayname',[refineType,' ','\eta_H'],'Marker','o','Markersize',10,'Linewidth',2,'Linestyle',linestyle)
% loglog(nrDoF,term5,'Displayname',[refineType,' ','\eta_{(1)}'],'Marker','<','Markersize',10,'Linewidth',2,'Linestyle',linestyle)
% loglog(nrDoF,term6,'Displayname',[refineType,' ','\eta_{(2)}'],'Marker','>','Markersize',10,'Linewidth',2,'Linestyle',linestyle)
set(gca,'Fontsize',18,'Xscale','log','Yscale','log')
legend('location','southwest')
