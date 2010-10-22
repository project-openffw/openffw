function p = ODRTestimate_PD(p)
% ||sigma-sigma_ell||_2^2 <= E(u_ell) + int W*(tau) dx
% with div tau = - f = - 1
%author: Lena Noack

%% INPUT
% length4ed = p.level(end).enum.length4ed;
lvl = size(p.level,2);
n4e = p.level(lvl).geom.n4e;
% n4ed = p.level(end).enum.n4ed;
% ed4e = p.level(end).enum.ed4e;
curLvl = length(p.level);
degree = 10;
nrNodes = p.level(end).nrNodes;
length4ed = p.level(end).enum.length4ed;
ed4e = p.level(end).enum.ed4e;

h_T = max(length4ed(ed4e)')';
%area4T = 0.25*h_T.^2; %works for triangles with 1 angle of 90° -> still to change;

% primal term: E(u_\ell)
intEnergy = integrate(n4e,lvl,degree,@getEnergy,p);
%EnMid_h = sum(intEnergy);
%fprintf('\nEnergy = %.15g\n',EnMid_h)

% dual term
W_Conj = integrate(n4e,lvl,degree,@getWDual,p);

nu = intEnergy + W_Conj;

%% OUTPUT
p.level(end).etaT = sqrt(abs(nu));
p.level(end).estimatedError = (norm(nu)).^0.5; %est. for ||sigma-sigma_h||

%% supply the discrete energy
function val = getEnergy(x,y,curElem,lvl,p)

energy_h = p.statics.energy_h;
val = energy_h(x,y,curElem,lvl,p);


function val = getWDual(x,y,curElem,lvl,p)

%q = -[x./2 y./2];               % div q = -1
sigma_h = p.statics.sigma_h;
q = sigma_h(x,y,curElem,lvl,p);
qAbs = ( q(:,1).^2 + q(:,2).^2 ).^(1/2);

WDual = p.problem.conjNonLinearFunc;
evalWDual = WDual(qAbs,curElem,lvl,p);

val(1,1,:) = evalWDual;
