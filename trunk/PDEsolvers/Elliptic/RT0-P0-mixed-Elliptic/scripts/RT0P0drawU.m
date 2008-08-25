function p = RT0P0drawU(p,lvl)
% draw displacement

% Copyright 2007 Jan Reininghaus, David Guenther
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


%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end

% set graphic options from structure p
drawInfo = loadField('p.params.output','drawInfo',p,true);
drawWalls = loadField('p.params.output','drawWalls',p,true);

% load geometry
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

% load discrete solution
u = p.level(lvl).u;

% load enumerated data
nrElems = p.level(lvl).nrElems;
nrDoF = p.level(lvl).nrDoF;

%% drawU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cla;
hold on

nrNodes4e = size(n4e,2);
Z = repmat(u,[nrNodes4e 1]);
coordX = reshape(c4n(n4e',1),nrNodes4e,nrElems);
coordY = reshape(c4n(n4e',2),nrNodes4e,nrElems);
patch(coordX,coordY,Z,Z);

X = zeros(4,nrNodes4e*nrElems);
Y = zeros(4,nrNodes4e*nrElems);
Z = zeros(4,nrNodes4e*nrElems);

if(drawWalls)
    curElem = 1:nrElems;
    curNodes = n4e(curElem,:);
    one = ones(size(curElem));
    curU = u(curElem);
    zero = zeros(size(curU));  
    
    coordX = reshape(c4n(curNodes(:,[1 2]),1),[nrElems 2]);
    coordY = reshape(c4n(curNodes(:,[1 2]),2),[nrElems 2]);  
    Z(:,nrNodes4e*(curElem-one)+one) = [zero;zero;curU;curU];
    X(:,nrNodes4e*(curElem-one)+one) = [coordX,coordX(:,[2 1])]';
    Y(:,nrNodes4e*(curElem-one)+one) = [coordY,coordY(:,[2 1])]';

    coordX = reshape(c4n(curNodes(:,[2 3]),1),[nrElems 2]);
    coordY = reshape(c4n(curNodes(:,[2 3]),2),[nrElems 2]);
    curU = u(curElem);
    Z(:,nrNodes4e*(curElem-one)+2*one) = [zero;zero;curU;curU];
    X(:,nrNodes4e*(curElem-one)+2*one) = [coordX,coordX(:,[2 1])]';
    Y(:,nrNodes4e*(curElem-one)+2*one) = [coordY,coordY(:,[2 1])]';

    coordX = reshape(c4n(curNodes(:,[3 1]),1),[nrElems 2]);
    coordY = reshape(c4n(curNodes(:,[3 1]),2),[nrElems 2]);
    curU = u(curElem);
    Z(:,nrNodes4e*(curElem-one)+3*one) = [zero;zero;curU;curU];
    X(:,nrNodes4e*(curElem-one)+3*one) = [coordX,coordX(:,[2 1])]';
    Y(:,nrNodes4e*(curElem-one)+3*one) = [coordY,coordY(:,[2 1])]';
end

patch(X,Y,Z,Z)
hold off
view(30,30)

if(drawInfo)
	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
	title('Discrete P0 Solution');
	grid on;
end