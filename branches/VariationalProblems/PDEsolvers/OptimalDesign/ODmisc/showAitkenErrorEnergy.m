function [p,errorEnergy] = showAitkenErrorEnergy(p,k)

energy_h = p.statics.energy_h;
nrLevels = size(p.level,2)

energy = zeros(nrLevels-1,1);
dof = zeros(nrLevels-1,1);


if k== 1
% OD Square uniform
eAitken = -0.015381759013978;
elseif k == 2
% OD Lshape uniform 
eAitken = -0.096307207969219;
elseif k == 3
% OD Octagon uniform
eAitken = -0.136825846977538;
elseif k == 4
% OD Slit uniform
eAitken = -0.146387380913813;
end
% 
for lvl = 2:nrLevels
    lvl
    nrElems = p.level(lvl).nrElems;
    dof(lvl-1) = p.level(lvl).nrDoF;
%     energy(lvl-1) = sum(integrate(p.level(lvl).geom.n4e,lvl,10,energy_h,p));
%     p.level(lvl).energy = energy(lvl-1);
    energy(lvl-1) = p.level(lvl).energy;
end

origEnergy = energy;
aitken = zeros(length(energy),1);

% while length(aitken)>2
%     aitken = getAitken(energy)
%     energy = aitken(end);
% end

% aitken
% aitken = aitken(1)

% limit = aitken(end-2)
% limit = eAitken;
% limit = -0.0149;
errorEnergy = abs(origEnergy-limit)

dof = dof(1:end-2);
errorEnergy = errorEnergy(1:end-2);

loglog(dof,errorEnergy,'kx-');

minDoF = 50;

I = find(dof >= minDoF);
nrDoF4lvl = dof(I);
err4lvl = errorEnergy(I);

nrDoF4lvl = log(nrDoF4lvl);
err4lvl = log(err4lvl);
approx = polyfit(nrDoF4lvl,err4lvl,1);

y = polyval(approx,nrDoF4lvl);

rate = -(y(end) - y(1)) / (nrDoF4lvl(end) - nrDoF4lvl(1))


function val = getAitken(x)

val = zeros(length(x)-2,1);

for k = 1:length(x)-2
    val(k) = x(k) - (x(k+1)-x(k))^2/...
                    (x(k+2)-2*x(k+1)+x(k));
end

