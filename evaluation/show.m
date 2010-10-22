function p = show(func,p)

mark = p.params.mark;
problem = p.params.problem.name;
pdeSolver = p.params.pdeSolver;
maxNrDoF = p.params.maxNrDoF;
nrLevels = p.level(end).level;
refineFirstLevel = p.params.modules.mark.refineFirstLevel;

showIteratively = loadField('p.params.output','showIteratively',p,false);
saveFigures = loadField('p.params.output','saveFigures',p,false);
format = loadField('p.params.output','format',p,'fig');


[token, rem] = strtok(func,'_');
drawFunc = str2func(token);


disp([func '_' pdeSolver '_' problem])

if(showIteratively)
	lvlVector = refineFirstLevel+1:nrLevels;
else
	lvlVector = nrLevels;
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
	pause(.5);

    
    % save files
	if( (saveFigures && ~strcmp(token,'drawError'))...
		|| (saveFigures && strcmp(token,'drawError') && lvl == nrLevels))
		warning off;
		sep = '_';
		mkdir('results');
		cd('results');
		directory = [pdeSolver,sep,problem,sep,num2str(maxNrDoF),sep,mark];
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
