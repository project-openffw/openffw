function p = showAitkenErrorEnergy(p)

energy_h = p.statics.energy_h;
nrLevels = size(p.level,2);

energy = zeros(nrLevels-1,1);
dof = zeros(nrLevels-1,1);

for lvl = 2:nrLevels
    nrElems = p.level(lvl).nrElems;
    dof(lvl-1) = p.level(lvl).nrDoF;
    n4e = p.level(lvl).geom.n4e;
%     intVal = integrate(n4e,lvl,19,energy_h,p);
%     energy(lvl-1) = sum(intVal)
%     p.level(lvl).energy = sum(intVal);
    energy(lvl-1) = p.level(lvl).energy;

end

origEnergy = energy
% aitken = zeros(length(energy),1);

% while length(aitken)>2
%     aitken = getAitken(energy);
%     energy = aitken;
% end

% aitken = -0.160477033993373;
aitken = -0.156810887275178;

aitken = aitken(end);

errorEnergy = abs(origEnergy-aitken(end))

% dof = dof(1:end-2);
% errorEnergy = errorEnergy(1:end-2);

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

val