function p = TWP2estimate(p)
% Lena Noack
%
%% Input
length4ed = p.level(end).enum.length4ed;
area4e = p.level(end).enum.area4e;
n4ed = p.level(end).enum.n4ed;
n4e = p.level(end).geom.n4e;
ed4e = p.level(end).enum.ed4e;
nrNodes = p.level(end).nrNodes;
degree = p.params.rhsIntegtrateExactDegree;
%p4n
%n4e
curLvl = size(p.level,2);

%% compute jump error and oscillations
eta4ed = length4ed.*integrate(n4ed,curLvl,2,@jumpError,p);
delta4ed = length4ed.*integrate(n4ed,curLvl,2,@jumpError,p);
[numbelems,dummy]=size(area4e);
for i=1:numbelems
    [k]=ed4e(i,:);
    delta4e(i,1) = delta4ed(k(1))+delta4ed(k(2))+delta4ed(k(3));
end

    osc = area4e.*integrate(n4e,curLvl,degree,@osc4_T,p);

%% compose error terms
estimatedError = sqrt( sum(eta4ed) + sum(osc) );

%% Output
%p.level(end).etaT  = eta4ed; 
p.level(end).deltaT  = sqrt(delta4e)+ area4e.*integrate(n4e,curLvl,degree,@f2,p); 
p.level(end).etaEd  = sqrt(eta4ed);
p.level(end).etaOsc = osc;
p.level(end).estimatedError = estimatedError;

return

%% supply jump error
function val = jumpError(x,y,curEdge,lvl,p)
sigma_h = p.statics.sigma_h;
e4ed = p.level(lvl).enum.e4ed;
normals4ed = p.level(lvl).enum.normals4ed;

elems = e4ed(curEdge,:);
normal = normals4ed(curEdge,:);

% assume kappa is the identity
kappa = repmat([1 0;0 1],[1 1 length(x)]);

evalSigma1 = sigma_h(x,y,elems(1),lvl,p)'*normal';
if elems(2) ~= 0
    % E is an interior edge
    evalSigma2 = sigma_h(x,y,elems(2),lvl,p)'*normal';
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
evalF = f(x,y,curElem,lvl,p);

val = zeros(1,1,length(x));
val(1,:,:) = (evalF - intMeanF).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = f2(x,y,curElem,lvl,p)

f = p.problem.f;
evalF = f(x,y,curElem,lvl,p);

val = zeros(1,1,length(x));
val(1,:,:) = (evalF).^2;

%% supply right hand side
function evalF = RHS(x,y,curElem,lvl,p)

f = p.problem.f;

evalF = f(x,y,curElem,lvl,p);
evalF = reshape(evalF,[1 1 length(x)]);