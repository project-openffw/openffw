function p = ODRTestimate(p)

% author: Lena Noack

%% INPUT
length4ed = p.level(end).enum.length4ed;
n4e = p.level(end).geom.n4e;
n4ed = p.level(end).enum.n4ed;
ed4e = p.level(end).enum.ed4e;
curLvl = length(p.level);
degree = loadField('p.params','rhsIntegrateExactDegree',p,5);

%% ESTIMATE
h_T = max(length4ed(ed4e)')';

%nu_T = h_T.^2.*integrate(n4e,curLvl,degree,@funcHandleResiduum,p);
nu_E = length4ed.*integrate(n4ed,curLvl,2,@jumpError,p);

%nu = sqrt( nu_T + 1/2*sum(nu_E(ed4e),2) );
nu = sqrt( sum(nu_E(ed4e),2) );

%% OUTPUT
p.level(end).etaT = nu;
% p.level(end).etaOsc = nu;
% p.level(end).etaEd = 0;
p.level(end).estimatedError = norm(nu);


function val = funcHandleResiduum(x,y,curElem,curLvl,p)
f = p.problem.f;
curf     = f(x,y,curElem,curLvl,p);

residuum = - curf(:);

%% RETURN
val(1,:,:) = (residuum.^2)';



function val = jumpError(x,y,curEdge,lvl,p)

sigma_h = p.statics.sigma_h;
e4ed = p.level(lvl).enum.e4ed;
normals4ed = p.level(lvl).enum.normals4ed;

elems = e4ed(curEdge,:);
normal = normals4ed(curEdge,:);

evalSigma1 = sigma_h(x,y,elems(1),lvl,p)*normal';

if elems(2) ~= 0
    % E is an interior edge
    evalSigma2 = sigma_h(x,y,elems(2),lvl,p)*normal';
else
    % E is a Dirichlet Edge
    evalSigma2 = evalSigma1;
end

val = zeros(1,1,length(x));
val(1,1,:) = sum((evalSigma1 - evalSigma2).^2,2);
