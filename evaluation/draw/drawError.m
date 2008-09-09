function p = drawError(p,lvl,errType)
% draws the convergence history.
%
% errType = estimatedError
%           H1semiError
%           L2error

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
lineStyle = loadField('p.params.output','lineStyle',p,'-');
myColor = loadField('p.params.output','myColor',p,'k');
marker = loadField('p.params.output','marker',p,'x');
minDoF = loadField('p.params.output','minDoF',p,ceil(p.level(end).nrDoF/17));
drawInfo = loadField('p.params.output','drawInfo',p,false);
plotGrid = loadField('p.params.output','plotGrid',p,true);
holdIt = loadField('p.params.output','holdIt',p,true);
name = loadField('p.params.output','name',p,[]);
fontSize = loadField('p.params.output','fontSize',p,[]);
setScales = loadField('p.params.output','setScales',p,false);
getConvergenceRate = loadField('p.params.output','getConvergenceRate',p,true);
drawConvergenceRate = loadField('p.params.output','drawConvergenceRate',p,false);
pointsForConvergenceRate = loadField('p.params.output','pointsForConvergenceRate',p,-1);
degree = loadField('p.params','errorIntegtrateExactDegree',p,19);
intMode = loadField('p.params','IntegrationMode',p,'elementwise');


%% draw routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for curLvl = 1:lvl
	eval(['dummy = p.level(curLvl).' errType ';'],'dummy = [];');
    if(isempty(dummy))
        if isempty(loadField('p.problem','u_exact',p))
            error('Cannot compute exact error because no exact solution was specified');
        end
        
        fprintf('\t Integrating, Level %2i of %2i, nrGaussPoints= %i',curLvl,p.level(end).level...
            , ceil((degree+1)/2 )^2 );
        time = cputime;
        if strcmpi(intMode,'elementwise')
            eval(['integrandFunc = p.statics.',errType,';']);
            err4e = integrate( p.level(curLvl).geom.n4e, curLvl,  degree, integrandFunc, p );
        else
            eval(['integrandFunc = p.statics.',errType,'Vectorised;']);
            err4e = integrateVectorised( p.level(curLvl).geom.n4e, curLvl,  degree, integrandFunc, p );
        end
        time = cputime - time;
        fprintf(', elapsed time = %.1g sec \n',time);
        eval(['p.level(curLvl).',errType,' = sqrt(sum(err4e));']);
        eval(['p.level(curLvl).',errType,'4e = sqrt(err4e);']);
    end
	err4lvl(curLvl) = eval(['p.level(curLvl).' errType]);
	nrDoF4lvl(curLvl) = p.level(curLvl).nrDoF;
end
	
I = find(nrDoF4lvl >= minDoF);
nrDoF4lvl = nrDoF4lvl(I);
err4lvl = err4lvl(I);

if(isempty(nrDoF4lvl))
	return
end

if(holdIt) 
    hold on; 
end
%if(isempty(name))
  if isfield(p.params,'legend') == 0;
      p.params.legend.name{1} = name;
  elseif isfield(p.params.legend,'name') == 0
      p.params.legend.name{1} = name;
  else
    p.params.legend.name{size(p.params.legend.name,2)+1} = name;
  end
%   p.params.colorIndex = p.params.colorIndex+1;
  p.params.colorIndex = 1+mod(p.params.colorIndex,size(p.params.colors,1)+1);
    loglog(nrDoF4lvl(1:end),err4lvl,p.params.colors(p.params.colorIndex,:));%,'Marker',marker,'linestyle',lineStyle);
    legend(p.params.legend.name)
    title(strcat('Convergence history - ',errType));
%    %loglog(nrDoF4lvl(1:end),err4lvl,'Marker',marker,'linestyle',lineStyle);
%else
%%     loglog(nrDoF4lvl(1:end),err4lvl,'Marker',marker,'linestyle',lineStyle,'DisplayName',name);
%    loglog(nrDoF4lvl(1:end),err4lvl);%,'Marker',marker,'linestyle',lineStyle);
%    %loglog(nrDoF4lvl(1:end),err4lvl,'Marker',marker,'linestyle',lineStyle);
%    legend('off')
%    legend('show')
%end
if(getConvergenceRate)
    dummy = log(nrDoF4lvl);
    approx = polyfit(dummy,log(err4lvl),1);
    y = polyval(approx,dummy);
    rate = -(y(end) - y(1)) / (dummy(end) - dummy(1));
    eval(['p.' errType '_rate = rate;']);
%    if(drawConvergenceRate)
%        loglog(exp(1).^dummy,exp(1).^y);
%    end
end
if(holdIt) 
    hold off; 
end

%set(gca,'XScale','log');
%set(gca,'YScale','log');

% set scale (produces problems when multiple graphs are drawn)
if(setScales)
    minXlog10 = floor(log10(min(nrDoF4lvl)));
    maxXlog10 = ceil(log10(max(nrDoF4lvl)));
    minYlog10 = floor(log10(min(err4lvl)));
    maxYlog10 = ceil(log10(max(err4lvl)));
    
    hold
    axis([10^minXlog10,10^maxXlog10,10^minYlog10,10^maxYlog10]);
%     hold off
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


