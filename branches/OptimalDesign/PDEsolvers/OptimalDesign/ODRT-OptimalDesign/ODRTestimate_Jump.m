function p = ODRTestimate_Jump(p)
%% INPUT
n4ed = p.level(end).enum.n4ed;
length4ed = p.level(end).enum.length4ed;

lvl = size(p.level,2);

degree = loadField('p.params','nonLinearExactIntegrateDegree',p,19);
%% estimation of the error
eta4ed = length4ed.*integrate(n4ed,lvl,degree,@intTangentJump,p);

%% compose the error terms
estimatedError = sqrt( sum(eta4ed) );

%% OUTPUT
p.level(end).etaEd = sqrt(eta4ed);
p.level(end).estimatedError = estimatedError;

%% supply tangent jump: int_E ([DW*(p)]*tau)^2
function val = intTangentJump(points,curEdge,lvl,p)

grad_h = p.statics.grad_h;
e4ed = p.level(lvl).enum.e4ed;
tangents4ed = p.level(lvl).enum.tangents4ed;

elems = e4ed(curEdge,:);
tangent = tangents4ed(curEdge,:)';

evalGrad_1 = grad_h(points,elems(1),lvl,p);
if elems(2) ~= 0
    evalGrad_2 = grad_h(points,elems(2),lvl,p);
else
    evalGrad_2 = evalGrad_1;
end

diffGrad = (evalGrad_1 - evalGrad_2)*tangent;
val = diffGrad.^2;
val = reshape(val,[1 1 length(points(:,1))]);
