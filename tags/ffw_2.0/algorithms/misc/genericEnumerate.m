function p = genericEnumerate(p)
% generate the generic data structures

% Copyright 2007 Jan Reininghaus, David Guenther, Andreas Byfut
%
% This file is part of FFW.
%
% FFW is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% FFW is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;
Db = p.level(end).geom.Db;
Nb = p.level(end).geom.Nb;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e4n = p.statics.enum.getE4n(n4e);
ed4n = p.statics.enum.getEd4n(e4n);
n4ed = p.statics.enum.getN4ed(ed4n);
ed4e = p.statics.enum.getEd4e(n4e,ed4n);
e4ed = p.statics.enum.getE4ed(n4ed,e4n);

DbEd = p.statics.enum.getDbEdges(Db,ed4n);
NbEd = p.statics.enum.getNbEdges(Nb,ed4n);

nrNodes = size(c4n,1);
nrElems = size(n4e,1);
nrEdges = size(n4ed,1);

area4e = p.statics.enum.getArea4e(c4n,n4e);
area4n = p.statics.enum.getArea4n(area4e,e4n);
length4ed = p.statics.enum.getLength4ed(c4n,n4ed);

midPoint4e = p.statics.enum.getMidPoint4e(c4n,n4e);
midPoint4ed = p.statics.enum.getMidPoint4ed(c4n,n4ed);

% outer-unit-normals of edges per element
normals4e = p.statics.enum.getNormals4e(c4n,n4e,length4ed,ed4e);

% unit-tangents of edges per element
tangents4e = p.statics.enum.getTangents4e(c4n,n4e,length4ed,ed4e);

% angles (in radians) per element
angles4e = p.statics.enum.getAngles4e(tangents4e);
angle4n = p.statics.enum.getAngle4n(angles4e,n4e,nrElems,nrNodes);

% outer-unit-normals of Neumann edges
normals4NbEd = p.statics.enum.getNormals4NbEd(c4n,Nb);

% outer-unit-normals of Dirichlet edges
normals4DbEd = p.statics.enum.getNormals4DbEd(c4n,Db);

% unit-normals per edges
normals4ed = p.statics.enum.getNormals4ed(normals4e,ed4e);

% unit-tangents per edges
tangents4ed = p.statics.enum.getNormals4ed(tangents4e,ed4e);

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
