function p =  drawErrorOnGrid(p,lvl,errType)
% draw error on the mesh

% Copyright 2007 Jan Reininghaus, David Guenther, 
%                Andreas Byfut, Joscha Gedicke
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

if(nargin < 3)
	error('You have to specify a the type of error ''(e.g. H1semiError)''');
end

if isempty(lvl)
	lvl = p.level(end).level;
end

% set graphic options from structure p
drawInfo = loadField('p.params.output','drawInfo',p,false);
drawWalls = loadField('p.params.output','drawWalls',p,true);

% load geometry
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

% load enumerated data
nrElems = p.level(lvl).nrElems;
nrDoF = p.level(lvl).nrDoF;
ed4e = p.level(end).enum.ed4e;

% integration parameters
degree = loadField('p.params','errorIntegtrateExactDegree',p,19);
intMode = loadField('p.params','IntegrationMode',p,'elementwise');


%% load or calulate desired error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch errType
    case {'H1semiError','L2error'}
        eval(['err4e = p.level(lvl).' errType '4e;'],'err4e = [];');
        if isempty(err4e)
            p = calcExactError(errType,lvl,p);
            err4e = eval(['p.level(lvl).',errType,'4e']);
        end
        
    case 'estimatedError'
        eval('err4e = p.level(lvl).etaT;','err4e = [];');
        if isempty(err4e)
            eval('err4Ed = p.level(lvl).etaEd;','err4Ed = [];');
            if ~isempty(err4Ed)
                etaEd = p.level(lvl).etaEd;
                err4e = zeros(nrElems,1);
                for curElem = 1:nrElems
                    err4e(curElem) = sum(etaEd(ed4e(curElem,:)));
                end
            else
                error('No estimated error available.');
            end;
        end
        
    otherwise
        string1 = ['err4e = p.level(lvl).',errType,'4e;'];
        string2 = 'error(''You have to specify a valid type of error \n (e.g. H1semiError, such that there is a field H1semiError4e in the structure.'');';        
        eval(string1,string2);                 
end


if (size(err4e,1) ~= nrElems) 
    u = err4e;
else
    u = err4e';
end


%% draw routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cla
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

    if nrNodes4e == 3
        coordX = reshape(c4n(curNodes(:,[3 1]),1),[nrElems 2]);
        coordY = reshape(c4n(curNodes(:,[3 1]),2),[nrElems 2]);
        curU = u(curElem);
        Z(:,nrNodes4e*(curElem-one)+3*one) = [zero;zero;curU;curU];
        X(:,nrNodes4e*(curElem-one)+3*one) = [coordX,coordX(:,[2 1])]';
        Y(:,nrNodes4e*(curElem-one)+3*one) = [coordY,coordY(:,[2 1])]';
    elseif nrNodes4e == 4
        coordX = reshape(c4n(curNodes(:,[3 4]),1),[nrElems 2]);
        coordY = reshape(c4n(curNodes(:,[3 4]),2),[nrElems 2]);
        curU = u(curElem);
        Z(:,nrNodes4e*(curElem-one)+3*one) = [zero;zero;curU;curU];
        X(:,nrNodes4e*(curElem-one)+3*one) = [coordX,coordX(:,[2 1])]';
        Y(:,nrNodes4e*(curElem-one)+3*one) = [coordY,coordY(:,[2 1])]';

        coordX = reshape(c4n(curNodes(:,[4 1]),1),[nrElems 2]);
        coordY = reshape(c4n(curNodes(:,[4 1]),2),[nrElems 2]);
        curU = u(curElem);
        Z(:,nrNodes4e*(curElem-one)+4*one) = [zero;zero;curU;curU];
        X(:,nrNodes4e*(curElem-one)+4*one) = [coordX,coordX(:,[2 1])]';
        Y(:,nrNodes4e*(curElem-one)+4*one) = [coordY,coordY(:,[2 1])]';
    end
end

patch(X,Y,Z,Z)
hold off

view(0,90)

if(drawInfo)
% 	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
% 	title('Discrete P0 Solution');
    title(sprintf([errType, ', l=%g, N=%g'],lvl,nrDoF));
	grid on;
end
