function p = show(func,p,curLvl)
% graphical output
% show is the main interface for graphical output in the ffw.
%
% input:
% func   - string that specifies the output
% p      - structure of the FFW
% curLvl - optional specifies which level to use, i.e. for the displacement
%          plot
%
% examples:  
% displacement plot: p = show('drawU',p);
%
% mesh plot        : p = show('drawGrid',p);
%
% convergence history for the estimated error:
%   p = show('drawError_estimatedError',p);
%
% convergence history for the energy error
%   p = show('drawError_H1semiError',p);
%
% convergence history for the L^2 error
%   p = show('drawError_L2error',p);

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


if(nargin < 3 || isempty(curLvl))
	curLvl = p.level(end).level;
end

if curLvl > p.level(end).level
    curLvl = p.level(end).level;
    warning('FFW:curLvl','curLvl is larger then number of levels!!!  Set to %d',curLvl);
end

%% INPUT
mark = p.params.mark;
problem = p.params.problem.name;
pdeSolver = p.params.pdeSolver;
maxNrDoF = p.params.maxNrDoF;
nrLevels = p.level(end).level;

showIteratively = loadField('p.params.output','showIteratively',p,false);
saveFigures = loadField('p.params.output','saveFigures',p,false);
format = loadField('p.params.output','format',p,'fig');
pauseTime = loadField('p.params.output','pauseTime',p,0.5);
%%

[token, rem] = strtok(func,'_');
drawFunc = str2func(token);


disp([func '_' pdeSolver '_' problem])

if(showIteratively)
	lvlVector = 1:nrLevels;
else
	lvlVector = curLvl;
end

for lvl = lvlVector
	if(~isempty(rem))
		[token2, rem2] = strtok(rem,'_');	
		if(isempty(rem2))
			p = drawFunc(p,lvl,token2);
		else
			[token3, rem3] = strtok(rem2,'_');
			if(isempty(rem3))
				p = drawFunc(p,lvl,token2,token3);
			else
				token4 = strtok(rem3,'_');
				p = drawFunc(p,lvl,token2,token3,token4);
			end
		end
	else
		p = drawFunc(p,lvl);
	end
	drawnow;

    if pauseTime >= 0
        pause(pauseTime);
    else
        pause;
    end

    
    % save files
	if( (saveFigures && ~strcmp(token,'drawError'))...
		|| (saveFigures && strcmp(token,'drawError') && lvl == nrLevels))
		warning off;
		sep = '_';
		mkdir('results');
		cd('results');
		directory = [pdeSolver,sep,problem,sep,num2str(maxNrDoF)];
		mkdir(directory);
		cd(directory)
		filename = [func sep directory sep mark sep 'Level' sep];
		filename = sprintf([filename '%03d'],lvl);
		saveas(gcf,filename,format);
		cd('..')
		cd('..')
		warning on;
	end
end
