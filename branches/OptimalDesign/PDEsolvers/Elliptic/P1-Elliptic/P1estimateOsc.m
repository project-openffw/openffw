function p = P1estimate(p)

%% Input
length4ed = p.level(end).enum.length4ed;
area4e = p.level(end).enum.area4e;
n4ed = p.level(end).enum.n4ed;
n4e = p.level(end).geom.n4e;
degree = p.params.rhsIntegtrateExactDegree;

curLvl = size(p.level,2);

getOsc4T = loadField('p.params','osc4T',p,false);
getOsc4EdPatch = loadField('p.params','osc4EdPatch',p,false);
getOsc4TPatch = loadField('p.params','osc4TPatch',p,false);

%% compute jump error and oscillations
eta4ed = length4ed.*integrate(n4ed,curLvl,2,@jumpError,p);

if getOsc4T
    osc = area4e.*integrate(n4e,curLvl,degree,@osc4_T,p);
elseif getOsc4EdPatch
    osc = length4ed.^2.*integrate(n4ed,curLvl,degree,@osc4_EdPatch,p);
elseif getOsc4TPatch
    osc = area4e.*integrate(n4e,curLvl,degree,@osc4_TPatch,p);
else
    osc = 0;
end

%% compose error terms
estimatedError = sqrt( sum(eta4ed) + sum(osc) );

%% Output
p.level(end).etaEd  = sqrt(eta4ed);
p.level(end).etaOsc = osc;
p.level(end).estimatedError = estimatedError;

%% supply jump error
function val = jumpError(x,y,curEdge,lvl,p)

g = p.problem.g;
sigma_h = p.statics.sigma_h;
e4ed = p.level(lvl).enum.e4ed;
NbEd = p.level(lvl).enum.NbEd;
normals4ed = p.level(lvl).enum.normals4ed;

elems = e4ed(curEdge,:);
normal = normals4ed(curEdge,:);

% assume kappa is the identity
kappa = repmat([1 0;0 1],[1 1 length(x)]);

evalSigma1 = sigma_h(x,y,elems(1),kappa,lvl,p)*normal';
if elems(2) ~= 0
    % E is an interior edge
    evalSigma2 = sigma_h(x,y,elems(2),kappa,lvl,p)*normal';
elseif ~isempty(find(curEdge == NbEd))
    % E is a Neumann edge
    evalSigma2 = g(x,y,(normal'*ones(1,length(x)))',p);
else
    % E is a Dirichlet Edge
    evalSigma2 = evalSigma1;
end

val = zeros(1,1,length(x));
val(1,1,:) = sum((evalSigma1 - evalSigma2).^2,2);

%% supply oscillations 4 triangle
function val = osc4_T(x,y,curElem,lvl,p)

f = p.problem.f;
degree = p.params.rhsIntegtrateExactDegree;
n4e = p.level(lvl).geom.n4e;
area4e = p.level(lvl).enum.area4e;

intMeanF = 1/area4e(curElem)*integrate(n4e(curElem,:),lvl,degree,@RHS,p);
evalF = f(x,y,p);

val = zeros(1,1,length(x));
val(1,:,:) = (evalF - intMeanF).^2;

%% supply oscillations 4 edge patch
function val = osc4_EdPatch(x,y,curEdge,lvl,p)

f = p.problem.f;
degree = p.params.rhsIntegtrateExactDegree;
n4e = p.level(lvl).geom.n4e;
e4ed = p.level(lvl).enum.e4ed;
area4e = p.level(lvl).enum.area4e;

% if E is a boundary edge, the corresponding boundary element is consider
% to be its own neighbour
elems = e4ed(curEdge,:);
if elems(2) == 0
    elems(2) = elems(1);
end

intF4E_1 = integrate(n4e(elems(1),:),lvl,degree,@RHS,p);
intF4E_2 = integrate(n4e(elems(2),:),lvl,degree,@RHS,p);
patchSize = area4e(elems(1)) + area4e(elems(2));

intMeanF = 1/patchSize*( intF4E_1+intF4E_2 );
evalF = f(x,y,p);

val = zeros(1,1,length(x));
val(1,:,:) = (evalF - intMeanF).^2;

%% supply oscillations 4 triangle patch
function val = osc4_TPatch(x,y,curElem,lvl,p)

f = p.problem.f;
degree = p.params.rhsIntegtrateExactDegree;
n4e = p.level(lvl).geom.n4e;
ed4e = p.level(lvl).enum.ed4e;
e4ed = p.level(lvl).enum.e4ed;
area4e = p.level(lvl).enum.area4e;

edges = ed4e(curElem,:);

% if T is a boundary element, then at least one edge is a boundary edge.
% The neighbour-index is then zero, which has to be removed.
neighbours = unique(e4ed(edges,:));
neighbours( neighbours == 0 ) = [];

intF4E = integrate(n4e(neighbours,:),lvl,degree,@RHS,p);
patchSize = sum(area4e(neighbours));

intMeanF = 1/patchSize*sum(intF4E);
evalF = f(x,y,p);

val = zeros(1,1,length(x));
val(1,:,:) = (evalF - intMeanF).^2;

%% supply right hand side
function evalF = RHS(x,y,curElem,lvl,p)

f = p.problem.f;

evalF = f(x,y,p);
evalF = reshape(evalF,[1 1 length(x)]);