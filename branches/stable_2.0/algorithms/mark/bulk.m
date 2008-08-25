function p = bulk(p)
% bulk criteria
% marks all edges that satisfy the bulk criteria

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

etaT = loadField('p.level(end)','etaT',p);
etaEd = loadField('p.level(end)','etaEd',p);
etaOsc = loadField('p.level(end)','etaOsc',p);

ed4e = p.level(end).enum.ed4e;
nrEdges = p.level(end).nrEdges;
nrElems = p.level(end).nrElems;

thetaEd = loadField('p.params.modules.mark.bulk','thetaEd',p,1/2);
thetaT = loadField('p.params.modules.mark.bulk','thetaT',p,1/2);
thetaOsc = loadField('p.params.modules.mark.bulk','thetaOsc',p,1/2);


%% Bulk Critera %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

markedEdges = false(nrEdges,1);
refineElems = false(nrElems,1);
refineElemsBisec5 = false(nrElems,1);

if (nnz(etaEd) ~= 0)
    [sortedEtaEd,I] = sort(etaEd,'descend');
    sumEtaEd = cumsum(sortedEtaEd.^2);
    k = find(sumEtaEd >= thetaEd * norm(etaEd,2)^2,1,'first');

    if(thetaEd == 1)
        markedEdges = true(nrEdges,1);
    else
        markedEdges(I(1:k)) = true;
    end
end

if (nnz(etaT) ~= 0)
    [sortedEtaT,I] = sort(etaT,'descend');
    sumEtaT = cumsum(sortedEtaT.^2);
    k = find(sumEtaT >= thetaT * norm(etaT,2)^2,1,'first');
    refineElems(I(1:k)) = true;
end

if (nnz(etaOsc) ~= 0)
    [sortedEtaOsc,I] = sort(etaOsc,'descend');
    sumEtaOsc = cumsum(sortedEtaOsc.^2);
    k = find(sumEtaOsc >= thetaOsc * norm(etaOsc,2)^2,1,'first');
    refineElemsBisec5(I(1:k)) = true;
end

markedEdges4e = ed4e((refineElems | refineElemsBisec5),:);
markedEdges(markedEdges4e(:)) = true;


%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.level(end).markedEdges = markedEdges';
p.level(end).refineElemsBisec5 = refineElemsBisec5;
p = p.statics.closure(p);
