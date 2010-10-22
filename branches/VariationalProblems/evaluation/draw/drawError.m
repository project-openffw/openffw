function p = drawError(p,lvl,errType)

if(nargin < 3)
	error('You have to specify the type of error ''(i.e. approxL2error)''');
end

if isempty(lvl)
	lvl = p.level(end).level;
end

maxNrDoF = p.level(end).nrDoF;

lineStyle = loadField('p.params.output','lineStyle',p,'-');
myColor = loadField('p.params.output','myColor',p,'k');
marker = loadField('p.params.output','marker',p,'x');
minDoF = loadField('p.params.output','minDoF',p,100);
drawInfo = loadField('p.params.output','drawInfo',p,false);
plotGrid = loadField('p.params.output','plotGrid',p,true);
holdIt = loadField('p.params.output','holdIt',p,true);
name = loadField('p.params.output','name',p,[]);
fontSize = loadField('p.params.output','fontSize',p,[]);
setScales = loadField('p.params.output','setScales',p,false);
getConvergenceRate = loadField('p.params.output','getConvergenceRate',p,true);
drawConvergenceRate = loadField('p.params.output','drawConvergenceRate',p,false);
degree = loadField('p.params','errorIntegtrateExactDegree',p,19);
refineFirstLevel = p.params.modules.mark.refineFirstLevel;

for curLvl = refineFirstLevel+1:lvl
	eval(['dummy = p.level(curLvl).' errType ';'],'dummy = [];');
    if(isempty(dummy))
        if isempty(loadField('p.problem','u_exact',p))
            error('Cannot compute exact error because no exact solution was specified');
        end
        eval(['integrandFunc = p.statics.',errType,';']);
        %fprintf('\t Integrating, Level %2i of %2i, nrGaussPoints= %i',curLvl,p.level(end).level...
        %    , ceil((degree+1)/2 )^2 );
        time = cputime;
        err4e = integrate( p.level(curLvl).geom.n4e, curLvl,  degree, integrandFunc, p );
        time = cputime - time;
        %fprintf(', elapsed time = %.1g sec \n',time);
        paramSquare = loadField('p.params','square',p,0); 
        probClass = loadField('p.problem','probClass',p,'LSquare'); %L^2 or L^(4/3)
        if strcmp(probClass,'TwoWell') % L^(4/3) instead of L^2
            if paramSquare==1
                eval(['p.level(curLvl).',errType,' = (sum(err4e)).^(3/2);']);
                eval(['p.level(curLvl).',errType,'4e = err4e.^(3/2);']);
            else
                eval(['p.level(curLvl).',errType,' = (sum(err4e)).^(3/4);']);
                eval(['p.level(curLvl).',errType,'4e = err4e.^(3/4);']);
            end
        else
            if paramSquare==1
                eval(['p.level(curLvl).',errType,' = sum(err4e);']);
                eval(['p.level(curLvl).',errType,'4e = err4e;']);
            else
                eval(['p.level(curLvl).',errType,' = sqrt(sum(err4e));']);
                eval(['p.level(curLvl).',errType,'4e = sqrt(err4e);']);
            end
        end
    end
	err4lvl(curLvl-refineFirstLevel) = eval(['p.level(curLvl).' errType]);
	nrDoF4lvl(curLvl-refineFirstLevel) = p.level(curLvl).nrDoF;
end
	
I = find(nrDoF4lvl > minDoF);
nrDoF4lvl = nrDoF4lvl(I);
err4lvl = err4lvl(I);

if(isempty(nrDoF4lvl))
	return
end

if(holdIt) hold all; end
if(isempty(name))    
    loglog(nrDoF4lvl(1:end),err4lvl,'Marker',marker,'linestyle',lineStyle);
else
    loglog(nrDoF4lvl(1:end),err4lvl,'Marker',marker,'linestyle',lineStyle,'DisplayName',name);
    legend('off')
    legend('show')
end
if(getConvergenceRate)
    dummy = log(nrDoF4lvl);
    approx = polyfit(dummy,log(err4lvl),1);
    y = polyval(approx,dummy);
    rate = -(y(end) - y(1)) / (dummy(end) - dummy(1));
    eval(['p.' errType '_rate = rate;']);
    if(drawConvergenceRate)
        loglog(exp(1).^dummy,exp(1).^y);
    end
end
if(holdIt) hold off; end

set(gca,'XScale','log');
set(gca,'YScale','log');

% set scale (produces problems when multiple graphs are drawn)
if(setScales)
    minXlog10 = floor(log10(min(nrDoF4lvl)));
    maxXlog10 = ceil(log10(max(nrDoF4lvl)));
    minYlog10 = floor(log10(min(err4lvl)));
    maxYlog10 = ceil(log10(max(err4lvl)));
    
    hold on
    axis([10^minXlog10,10^maxXlog10,10^minYlog10,10^maxYlog10]);
    hold off
end
if(~isempty(fontSize))
    set(gca,'FontSize',fontSize);
end

if(drawInfo)
    xlabel(sprintf('Nr of degrees of freedom'));
    title(errType)
end
if(plotGrid) 
	grid on
else
    grid off
end


