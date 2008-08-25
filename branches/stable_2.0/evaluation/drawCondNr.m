function p = drawCondNr(p,lvl)
% draws the condition number of the stiffness matrix

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


if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end

%% input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set graphic options from structure p
lineStyle = loadField('p.params.output','lineStyle',p,'-');
myColor = loadField('p.params.output','myColor',p,'k');
marker = loadField('p.params.output','marker',p,'x');
minDoF = loadField('p.params.output','minDoF',p,ceil(p.level(end).nrDoF/17));
drawInfo = loadField('p.params.output','drawInfo',p,false);
plotGrid = loadField('p.params.output','plotGrid',p,true);
name = loadField('p.params.output','name',p,[]);
holdIt = loadField('p.params.output','holdIt',p,true);


%% draw routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = conditionNumber(p);

nrDoF4lvl = zeros(lvl,1);
cond4lvl = zeros(lvl,1);
for curLvl = 1:lvl
    nrDoF4lvl(curLvl) = p.level(curLvl).nrDoF;
    cond4lvl(curLvl) = p.level(curLvl).conditionNr;
end

I = find(nrDoF4lvl > minDoF);
nrDoF4lvl = nrDoF4lvl(I);
cond4lvl = cond4lvl(I);

if ~isempty(nrDoF4lvl)

    if holdIt 
        hold all; 
    end

    if(isempty(name))
        loglog(nrDoF4lvl(1:end),cond4lvl,'Marker',marker,'linestyle',lineStyle,'Color',myColor);
    else
        loglog(nrDoF4lvl(1:end),cond4lvl,'Marker',marker,'linestyle',lineStyle,'Color',myColor,'DisplayName',name);
        legend('off')
        legend('show')
    end
    set(gca,'XScale','log');
    set(gca,'YScale','log');

    if(drawInfo)
        xlabel(sprintf('Nr of degrees of freedom'));
        title('Condition number of global energy matrix')
    end

    if plotGrid
        grid on
    else
        grid off
    end

    hold off
    
end


%% local function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = conditionNumber(p)

level = p.level(end).level;

for curLvl = 1:level
    freeNodes = p.level(curLvl).enum.freeNodes;
    A = p.level(curLvl).A;
    p.level(curLvl).conditionNr = condest(A(freeNodees,freeNodes));
end
    
    

