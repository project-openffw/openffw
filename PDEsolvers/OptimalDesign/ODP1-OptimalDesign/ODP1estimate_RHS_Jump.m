function p = ODP1estimate_RHS_Jump(p)
%author: Lena Noack

p.statics.Ju_h = @getJu_h;
%p.level(end).Juh = @getJuh;
p.statics.Juh = @getJuh;

%% INPUT 
n4e = p.level(end).geom.n4e;
n4ed = p.level(end).enum.n4ed;
e4ed = p.level(end).enum.e4ed;
curLvl = length(p.level);
degree = loadField('p.params','rhsIntegrateExactDegree',p,5);

jumpS4ed = integrate(n4ed,curLvl,degree,@integrand1,p);  %int_E [sigma_h]nu (u_h-Ju_h) ds         

%for curEdge=1:length(e4ed) %int_{\omega_E} (RHS(u_h) - RHS(Ju_h)) dx
%    if e4ed(curEdge,2) ~= 0
%       intRhs4ed(curEdge) = integrate(n4e(e4ed(curEdge,1),:),curLvl,degree,@integrand2,p) + ...
%                         integrate(n4e(e4ed(curEdge,2),:),curLvl,degree,@integrand2,p);
%    else    
%       intRhs4ed(curEdge) = integrate(n4e(e4ed(curEdge,1),:),curLvl,degree,@integrand2,p);
%    end
%end

%eta4ed = sqrt(jumpS4ed + intRhs4ed);
eta4ed = sqrt(jumpS4ed);

%% OUTPUT 
p.level(end).jump4ed = sqrt(jumpS4ed);
p.level(end).etaEd = eta4ed;
p.level(end).estimatedError = norm(eta4ed,2);


%Supply integrand: [sigma_h]nu (u_h-Ju_h)
function val=integrand1(x,y,curEdge,lvl,p)

sigma_h = p.statics.sigma_h;
e4ed = p.level(end).enum.e4ed;
length4ed = p.level(end).enum.length4ed;
midPoint4e = p.level(end).enum.midPoint4e;
nrElems = p.level(end).nrElems;

normals4ed = p.level(lvl).enum.normals4ed;
nrEdges = p.level(end).nrEdges;

%nrElems
%nrEdges
%size(curEdge)

u_h = p.statics.u_h;
Ju_h = p.statics.Ju_h;




%% Compute the jump error of the stress
evalSigma = zeros(nrElems,2);
for curElem = 1:nrElems
    evalSigma(curElem,:) = sigma_h(midPoint4e(curElem,1),midPoint4e(curElem,1),...
                                    curElem,lvl,p);
end


phU = evalSigma(:,1);
phV = evalSigma(:,2);

e4ed = sort(e4ed,2);
[I,dontUse] = find(e4ed == 0);
e4ed(I,1) = e4ed(I,2);

ph4edU1 = phU(e4ed(curEdge,1));
ph4edU2 = phU(e4ed(curEdge,2));
ph4edU = ph4edU1 - ph4edU2;

ph4edV1 = phV(e4ed(curEdge,1));
ph4edV2 = phV(e4ed(curEdge,2));
ph4edV = ph4edV1 - ph4edV2;

DiffU = u_h(midPoint4e(e4ed(curEdge,1),1),midPoint4e(e4ed(curEdge,1),1),e4ed(curEdge,1),lvl,p)+...
        u_h(midPoint4e(e4ed(curEdge,1),1),midPoint4e(e4ed(curEdge,1),1),e4ed(curEdge,2),lvl,p);
DiffJuh = Ju_h(midPoint4e(e4ed(curEdge,1),1),midPoint4e(e4ed(curEdge,1),1),e4ed(curEdge,1),lvl,p)+...
          Ju_h(midPoint4e(e4ed(curEdge,1),1),midPoint4e(e4ed(curEdge,1),1),e4ed(curEdge,2),lvl,p);


normal = zeros(nrEdges,2);
%for curEdge = 1:nrEdges
    normal(curEdge,:) = normals4ed(curEdge,:);
%end

%size(sqrt((ph4edU.*normal(curEdge,1)).^2 + (ph4edV.*normal(curEdge,2)).^2))
%size((DiffU - DiffJuh))
%val = zeros(1,1,length(x));
val(1,1,:) = sqrt((ph4edU.*normal(curEdge,1)).^2 + (ph4edV.*normal(curEdge,2)).^2).*... ; %[sigma_h]nu
             abs(DiffU - DiffJuh);



% Supply integrand: RHS(u_\ell) - RHS(Ju_\ell)
function val=integrand2(x,y,curElem,lvl,p)

val = zeros(length(x));
%Noch dazu: RHS(u_h)-RHS(Ju_h) für nichtlineare Seite -> wie???
%f = p.problem.f;
%evalF = f(x,y,curElem,lvl,p);
%evalF = reshape(evalF,[1 1 length(x)]);

%% supply weighted approximation operator Ju_h for current element
function Ju_h = getJu_h(x,y,curElem,lvl,p)

%Ju_h4n = p.level(lvl).Juh; 
Ju_h4n = p.statics.Juh;
n4e = p.level(lvl).geom.n4e;

%Ju_h = Ju_h4n(n4e(curElem,:));
Ju_h = sum(Ju_h4n(x,y,n4e(curElem,:),lvl,p));

%% supply weighted approximation operator Ju_h
function Juh = getJuh(x,y,curNode,curLvl,p)

%curLvl = length(p.level);
n4e = p.level(curLvl).geom.n4e;
ed4e = p.level(curLvl).enum.ed4e;
e4ed = p.level(curLvl).enum.e4ed;
area4e = p.level(curLvl).enum.area4e;
midPoint4e = p.level(curLvl).enum.midPoint4e;

%curNode
for k=1:size(curNode,2)
  neighbours=neighbour(curNode(1,k),curLvl,p);
  if (size(neighbours,2)>0)
    basisU = p.statics.basisU;
%    neighbours
%    size(n4e)
    phiE = basisU(x,y,neighbours,curLvl,p);
    uh = p.statics.u_h;
    evalJuh = 0;
    for nb = 1:length(neighbours)
      uEs(nb) = uh(midPoint4e(nb,1),midPoint4e(nb,1),neighbours(nb),curLvl,p);
      phiEs(nb) = abs(phiE(nb));
      evalJuh = evalJuh + area4e(nb)*uEs(nb)*phiEs(nb);
    end
%size(area4e(neighbours))
%size(uEs')
%size(phiEs')
    Juh_t(k) = evalJuh;%sum(area4e(neighbours).*uEs'.*phiEs');
  else
    Juh_t(k) = 0;
  end
end
patchEs=patch(curNode,curLvl,p);
basisU = p.statics.basisU;

Iz4n = p.level(curLvl).enum.Iz4n;

psiE=[];
for k=1:size(curNode,2)
  Iz=[];
  [dontUse,Iz]=find(Iz4n==curNode(1,k));
  if (length(Iz')>=1)
    for node = 1:length(Iz')
      if (length(patch(node,curLvl,p))>=1)
        psiE = psiE + basisU(x,y,patch(node,curLvl,p),curLvl,p)';  
      else psiE=0;
      end
    end
  else psiE=0;
  end
  psiEs = abs(psiE);
  Juh_b(k) = sum(area4e(patch(curNode(1,k),curLvl,p))*psiEs);
  if (Juh_b(k) == 0)
    Juh_b(k) = 1;
  end
end

%size(Juh_t)
%size(Juh_b)
Juh = Juh_t/Juh_b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = neighbour(curNode,lvl,p)

p4n = p.level(lvl).enum.p4n;
e4n = p.level(lvl).enum.e4n;

[I,neighbours]=find(p4n(curNode,:));
neighbours( neighbours == 0 ) = [];

val = neighbours;

function val = patch(curNode,lvl,p)

p4n = p.level(lvl).enum.wp4n;

[I,patchE]=find(p4n(curNode,:));
patchE( patchE == 0 ) = [];

val = patchE;

