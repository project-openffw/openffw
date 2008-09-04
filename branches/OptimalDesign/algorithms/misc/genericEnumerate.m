function p = genericEnumerate(p)

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;
Db = p.level(end).geom.Db;
Nb = p.level(end).geom.Nb;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e4n = getE4n(n4e);
ed4n = getEd4n(e4n);
n4ed = getN4ed(ed4n);
ed4e = getEd4e(n4e,ed4n);
e4ed = getE4ed(n4ed,e4n);

DbEd = getDbEdges(Db,ed4n);
NbEd = getNbEdges(Nb,ed4n);

nrNodes = size(c4n,1);
nrElems = size(n4e,1);
nrEdges = size(n4ed,1);

area4e = getArea4e(c4n,n4e);
area4n = getArea4n(area4e,e4n);
length4ed = getLength4ed(c4n,n4ed);

midPoint4e = getMidPoint4e(c4n,n4e);
midPoint4ed = getMidPoint4ed(c4n,n4ed);

% outer-unit-normals of edges per element
normals4e = getNormals4e(c4n,n4e,length4ed,ed4e);

% unit-tangents of edges per element
tangents4e = getTangents4e(c4n,n4e,length4ed,ed4e);

% angles (in radians) per element
angles4e = getAngles4e(tangents4e);
angle4n = getAngle4n(angles4e,n4e,nrElems,nrNodes);

% outer-unit-normals of Neumann edges
normals4NbEd = getNormals4NbEd(c4n,Nb);

% outer-unit-normals of Dirichlet edges
normals4DbEd = getNormals4DbEd(c4n,Db);

% unit-normals per edges
normals4ed = getNormals4ed(normals4e,ed4e);

% unit-tangents per edges
tangents4ed = getNormals4ed(tangents4e,ed4e);

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).enum.ed4n = ed4n;
p.level(end).enum.ed4e = ed4e;
p.level(end).enum.e4n = e4n;
p.level(end).enum.n4ed = n4ed;
p.level(end).enum.e4ed = e4ed;
p.level(end).enum.DbEd = DbEd;
p.level(end).enum.NbEd = NbEd;
p.level(end).enum.area4e = area4e;
p.level(end).enum.area4n = area4n;
p.level(end).enum.length4ed = length4ed;
p.level(end).enum.normals4e = normals4e;
p.level(end).enum.normals4ed = normals4ed;
p.level(end).enum.tangents4e = tangents4e;
p.level(end).enum.tangents4ed = tangents4ed;
p.level(end).enum.angles4e = angles4e;
p.level(end).enum.angle4n = angle4n;
p.level(end).enum.normals4NbEd = normals4NbEd;
p.level(end).enum.normals4DbEd = normals4DbEd;
p.level(end).enum.midPoint4e = midPoint4e;
p.level(end).enum.midPoint4ed = midPoint4ed;
p.level(end).nrNodes = nrNodes;
p.level(end).nrElems = nrElems;
p.level(end).nrEdges = nrEdges;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
