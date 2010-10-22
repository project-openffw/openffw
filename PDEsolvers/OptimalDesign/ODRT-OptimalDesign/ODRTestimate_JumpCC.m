function p = ODRTestimate_JumpCC(p)
%author: Lena Noack
% ( h_T^2 ||f +div sigma_h||^2_(L^2) + 0.5 sum h_E ||[sigma_h]nu||_{L^2} )^{1/2}

%% INPUT
length4ed = p.level(end).enum.length4ed;
n4e = p.level(end).geom.n4e;
n4ed = p.level(end).enum.n4ed;
ed4e = p.level(end).enum.ed4e;
curLvl = length(p.level);
degree = loadField('p.params','rhsIntegrateExactDegree',p,5);

grad_h = p.statics.GradU_h;  % grad_h=p_h
e4ed = p.level(end).enum.e4ed;
length4ed = p.level(end).enum.length4ed;
midPoint4e = p.level(end).enum.midPoint4e;
nrElems = p.level(end).nrElems;
lvl = length(p.level);

normals4ed = p.level(lvl).enum.normals4ed;
nrEdges = p.level(end).nrEdges;

h_T = max(length4ed(ed4e),[],2);

mu_T = integrate(n4e,curLvl,2*degree,@Residuum,p);
mu_E = integrate(n4ed,curLvl,degree,@integrand,p);
mu = sqrt( h_T.^2.*mu_T + 1/2*h_T.*sum(mu_E(ed4e),2) );

%% OUTPUT 
p.level(end).etaT = mu;
p.level(end).estimatedError = sqrt(norm(mu));  %est. for ||sigma-sigma_h||
%p.level(end).estimatedError = norm(mu);  %est. for ||sigma-sigma_h||^2

function val = Residuum(x,y,curElem,curLvl,p)
f = p.problem.f;
curf     = f(x,y,curElem,curLvl,p);

D2W = p.problem.nonLinearExactSecDer;
grad_h = p.statics.GradU_h;
evalGrad = grad_h(x,y,curElem,curLvl,p);
absGrad = ( evalGrad(:,1).^2 + evalGrad(:,2).^2 ).^(1/2);
evalD2W = D2W(absGrad,curElem,curLvl,p); %sigma_h = D2W

% still to change: full derivative D2W -> D2W+DW (see prob. def.)

residuum = - curf(:);

sumRes = residuum + evalD2W;

val(1,:,:) = (sumRes.^2)';

function val = integrand(x,y,curEdge,lvl,p)

%sigma_h = p.statics.sigma_h; %sigma_h = DW(p_h)
%e4ed = p.level(lvl).enum.e4ed;
%normals4ed = p.level(lvl).enum.normals4ed;

%elems = e4ed(curEdge,:);
%normal = normals4ed(curEdge,:);

%evalSigma1 = sigma_h(x,y,elems(1),lvl,p)*normal'; %sigma_h|T_1 * nu_E

%if elems(2) ~= 0
    % E is an interior edge
%    evalSigma2 = sigma_h(x,y,elems(2),lvl,p)*normal'; %sigma_h|T_2 * nu_E
%else
    % E is a Dirichlet Edge
%    evalSigma2 = evalSigma1;
%end

%val = zeros(1,1,length(x));
%val(1,1,:) = sum((evalSigma1 - evalSigma2).^2,2); %sum( (sigma_h|T_1 - sigma_h|T_2 ) * nu_E )^2

%sigma_h = p.statics.sigma_h;
grad_h = p.statics.grad_h;
DW = p.problem.nonLinearExactDer;                       %DW(Du)/|Du|

e4ed = p.level(lvl).enum.e4ed;
normals4ed = p.level(lvl).enum.normals4ed;

elems = e4ed(curEdge,:);
normal = normals4ed(curEdge,:);

%evalSigma1 = sigma_h(x,y,elems(1),lvl,p)*normal';
evalGrad = grad_h(x,y,elems(1),lvl,p);
absGrad = ( evalGrad(:,1).^2 + evalGrad(:,2).^2 ).^(1/2);
evalSigma1 = DW(absGrad,elems(1),lvl,p).*absGrad;

if elems(2) ~= 0
    % E is an interior edge
%    evalSigma2 = sigma_h(x,y,elems(2),lvl,p)*normal';
    evalGrad = grad_h(x,y,elems(2),lvl,p);
    absGrad = ( evalGrad(:,1).^2 + evalGrad(:,2).^2 ).^(1/2);
    evalSigma2 = DW(absGrad,elems(2),lvl,p).*absGrad;
else
    % E is a Dirichlet Edge
    evalSigma2 = evalSigma1;
end

val = zeros(1,1,length(x));
val(1,1,:) = sum((evalSigma1 - evalSigma2).^2,2);

