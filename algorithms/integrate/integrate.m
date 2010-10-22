function val = integrate( parts, curLvl,  degree, integrand, p )
%INTEGRATE Integrates function handle integrand over given domain 1D or 2D
% integrand must have the form integrand(x,y,curPart,curLvl,p) and the
% return value must have the dimension [ n m nrGaussPoints ]

%% INPUT

area4e = p.level(curLvl).enum.area4e;
length4ed = p.level(curLvl).enum.length4ed;
e4n = p.level(curLvl).enum.e4n;
ed4n = p.level(curLvl).enum.ed4n;
nrParts = size(parts,1);


%% Number of Gauss Points

nrGaussPoints = ceil((degree+1)/2);

%% 1D or 2D
if size(parts,2) == 2
    % 1D
    [refGaussPoints,refGaussWeights] = getGaussPoints(nrGaussPoints);
    indices = rowaddr(ed4n,parts(:,1),parts(:,2));
    weightFactor = length4ed(indices);

elseif size(parts,2) == 3
    % 2D triangle
    [refGaussPoints,refGaussWeights] = getConProdGaussPoints(nrGaussPoints);
    indices = rowaddr(e4n,parts(:,1),parts(:,2));
    weightFactor = 2*area4e(indices);
    nrGaussPoints = nrGaussPoints^2;

elseif size(parts,2) == 4
    % 2D square
    [refGaussPoints,refGaussWeights] = getCompositeGaussPoints(nrGaussPoints);
    indices = rowaddr(e4n,parts(:,1),parts(:,2));
    weightFactor = area4e(indices);
    nrGaussPoints = nrGaussPoints^2;

else
    error('Could not integrate because no Domain is specified.');
end

%% Integrate
for curPart = 1:nrParts

    % transformation
    curCoords = ref2arbitrary(parts(curPart,:),refGaussPoints,curLvl,p);
    % evaluate function handle eval must be [n m nrGaussPoints]
    integrandVal = integrand(curCoords(:,1), curCoords(:,2), indices(curPart), curLvl, p);
    if(curPart == 1)
        dim = size(integrandVal);
        val = zeros(nrParts,dim(1),dim(2));
    end
    % weighted sum
    curWeights = weightFactor(curPart)*refGaussWeights;
    curWeights = reshape( (curWeights*ones(1,dim(1)*dim(2)))',[dim(1),dim(2),nrGaussPoints]);
    val(curPart,:,:) =  sum(curWeights.*integrandVal,3);

end

%% ref2arbitrary
function val = ref2arbitrary(part,refGaussPoints,curLvl,p)

c4n = p.level(curLvl).geom.c4n;
coords = c4n(part,:);
nrGaussPoints = size(refGaussPoints,1);

%1D or 2D
if length(part) == 2
    %1D ref = conv(0,1)
    c1 = (coords(1,:)'*ones(1,nrGaussPoints))';
    c2 = (coords(2,:)'*ones(1,nrGaussPoints))';

    val = c1 + [refGaussPoints(:,1),refGaussPoints(:,1)].*(c2-c1);

elseif length(part) == 3
    % 2D ref = conv{ (0,0),(1,0),(0,1)}
    c1 = (coords(1,:)'*ones(1,nrGaussPoints))';
    c2 = (coords(2,:)'*ones(1,nrGaussPoints))';
    c3 = (coords(3,:)'*ones(1,nrGaussPoints))';

    val = c1 + [refGaussPoints(:,1),refGaussPoints(:,1)].*(c2-c1) + ...
        [refGaussPoints(:,2),refGaussPoints(:,2)].*(c3-c1);

elseif length(part) == 4
    % 2D ref = conv{(0,0),(1,0),(0,1),(1,1)}
    c1 = (coords(1,:)'*ones(1,nrGaussPoints))';
    c2 = (coords(2,:)'*ones(1,nrGaussPoints))';
    c3 = (coords(3,:)'*ones(1,nrGaussPoints))';
    c4 = (coords(4,:)'*ones(1,nrGaussPoints))';

    val = c1 + [refGaussPoints(:,1),refGaussPoints(:,1)].*(c2-c1) + ...
        [refGaussPoints(:,2),refGaussPoints(:,2)].* ...
        (c4-c1 + [refGaussPoints(:,1),refGaussPoints(:,1)].*(c3-c4-c2+c1));
end
