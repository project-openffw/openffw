% draw convergence history
function evalExamples(p)
% Copyright 2007 David Guenther
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
%
%
%

% string = ['structA' 'Jump'];
% string = ['ODRTGOALstruct' 'Slit']
% load(string)
k = 4; 
% nrStructs = length(A.struct);

%square 0.027983242756785
% lshape 0.160674590782529
%octagon 0.257317403528310
%slit 0.247295558950406

% for k = 1:nrStructs
%    p = A.struct(k).p;
%    
   if k < 3
       limit = 0.027983242756785;
   elseif k < 5
       limit = 0.160674590782529;
   elseif k < 7
       limit = 0.257317403528310;
   elseif k < 9
       limit = 0.247295558950406;
   end

   for lvl = 2:length(p.level);
%        energy(lvl-1) = sum(integrate(p.level(lvl).geom.n4e,lvl,10,@getConjEnergy,p));
       nrDoF(lvl-1) = p.level(lvl).nrDoF;
   end
   
%    extrapolEnergy = getAitken(energy(4:end))'
   
%    if (mod(k,2)==1)
%     limit = extrapolEnergy(end)
%    end
   
%    errorEnergy = abs(energy - limit);
   
%    p = drawFigures(p,k);
%    hold all
%    loglog(nrDoF,errorEnergy,'o-');
%    clear energy nrDoF
% end

% for k = 1:nrStructs
%    p = A.struct(k).p;
   
%    if mod(k,2) == 1
%        figure
%        ODshowVolumeFraction(p,length(p.level),true);
%        colorbar
%        string = ['volumeFrac' num2str(k)];
%        saveas(gcf,string,'fig');
%        saveas(gcf,string,'epsc2');
%    end
   
%    if mod(k,2) == 0
%        nrLvl = length(p.level);
%        figure
%        drawGrid(p,1);
%        string = ['Eta2grid1Geom' num2str(k)];
%        set(gca,'FontSize',18);
%        saveas(gcf,string,'fig');
%        saveas(gcf,string,'epsc2');
%        
%        figure
%        drawGrid(p,ceil(nrLvl*1/4));
%        set(gca,'FontSize',18);
%        string = ['Eta2grid2Geom' num2str(k)];
%        saveas(gcf,string,'fig');
%        saveas(gcf,string,'epsc2');
%        
%        figure
%        drawGrid(p,ceil(nrLvl*2/4));
%        set(gca,'FontSize',18);
%        string = ['Eta2grid3Geom' num2str(k)];
%        saveas(gcf,string,'fig');
%        saveas(gcf,string,'epsc2');
%        
%        figure
%        drawGrid(p,ceil(nrLvl*3/4));
%        set(gca,'FontSize',18);
%        string = ['Eta2grid4Geom' num2str(k)];
%        saveas(gcf,string,'fig');
%        saveas(gcf,string,'epsc2');
%          dummy1 = [1 ceil(nrLvl*1/4) ceil(nrLvl*2/4) ceil(nrLvl*3/4)]
%          dummy2 = [p.level(1).nrDoF p.level(ceil(nrLvl*1/4)).nrDoF p.level(ceil(nrLvl*2/4)).nrDoF p.level(ceil(nrLvl*3/4)).nrDoF]
%    end
%    
% end

% for k = 1:2
%    p = A.struct(k).p;
   
    volumeFraction = p.statics.volumeFraction;
    q = p;
    for lvl = 2:length(p.level)
        c4n = p.level(lvl).geom.c4n;
        n4e = p.level(lvl).geom.n4e;
        coordsX = reshape(c4n(n4e,1),[],3)';
        coordsY = reshape(c4n(n4e,2),[],3)';
        nrElems = p.level(lvl).nrElems;
        val = zeros(3,nrElems);

        for curElem = 1:nrElems
            val(:,curElem) = volumeFraction(coordsX(:,curElem),coordsY(:,curElem),curElem,lvl,p);
        end

        dummy = sum(val,1)';
        I = find((dummy ~= 3) & (dummy ~= 0));
        area4e = p.level(lvl).enum.area4e;
        areaContactZone(lvl-1) = sum(area4e(I));
        averageError = integrate(n4e,lvl,10,@phminusAph,p);
        p.level(lvl).phminusAph = sqrt(sum(averageError));
        p.level(lvl).phminusAph4e = sqrt(averageError);
        errorAverage(lvl-1) = sqrt(sum(averageError));

%         I = find((dummy == 3));
%         I = find((dummy == 0));
%         ed4e = p.level(lvl).enum.ed4e;
%         edges = unique(ed4e(I,:));
%         
%         e4ed = p.level(lvl).enum.e4ed;
%         
%         elems = e4ed(edges,:);
%         elems(elems(:,2) == 0,2) = elems(elems(:,2) == 0,1);
%         
%         elems = reshape(dummy(elems),size(elems,1),2);
%         index1 = find(sum(elems,2) == 0);
%         index2 = find(sum(elems,2) == 6);
%         
%         etaEd = p.level(lvl).etaEd;
%         
%         innerError = etaEd(edges([index1;index2]));

%         etaProj = q.level(lvl).etaProj4e;
%         innerError = etaProj(I);
        
%         estimatedError = sqrt(sum(innerError.^2));
%         q.level(lvl).estimatedError = estimatedError;
    end
    areaContactZone
%     limitArea = getAitken(totalArea(1:end))
%     q = show('drawError_estimatedError',q);
%     clear p q;
% end
% %%
firstIndex = find(nrDoF > 500,1,'first');

index = firstIndex:length(p.level)-1;
nrDoF = nrDoF(index);

term2 = areaContactZone(index).^(1/2);
term3 = errorAverage(index);
term4 = ((areaContactZone(index) + errorAverage(index)).*errorAverage(index)).^(1/2);

hold all
% loglog(nrDoF,term1,'Displayname',[refineType,' ','(E(\sigma_{\theta,h})-E(\sigma))^{1/2}'],'Marker','s','Markersize',10,'Linewidth',2,'Linestyle',linestyle)
loglog(nrDoF,term2,'Displayname',['\xi=0.8',' ','Area Contact Zone'],'Marker','d','Markersize',10,'Linewidth',2,'Linestyle','--')
loglog(nrDoF,term3,'Displayname',['\xi=0.8',' ','||\sigma_{h,\theta}-A\sigma_{h,\theta}||_{L^2}'],'Marker','p','Markersize',10,'Linewidth',2,'Linestyle','--')
loglog(nrDoF,term4,'Displayname',['\xi=0.8',' ','\eta_H'],'Marker','o','Markersize',10,'Linewidth',2,'Linestyle','--')
% loglog(nrDoF,term5,'Displayname',[refineType,' ','\eta_{(1)}'],'Marker','<','Markersize',10,'Linewidth',2,'Linestyle',linestyle)
% loglog(nrDoF,term6,'Displayname',[refineType,' ','\eta_{(2)}'],'Marker','>','Markersize',10,'Linewidth',2,'Linestyle',linestyle)
set(gca,'Fontsize',18,'Xscale','log','Yscale','log')
legend('location','southwest')

%% supply the error ||p_{h,\epsilon)-Ap_{h,\epsilon)||_L^2
function val = phminusAph(x,y,curElem,lvl,p)

sigma_h = p.statics.sigma_h;
Asigma_h = p.statics.Ap_h;

Asigmah = Asigma_h(x,y,curElem,lvl,p);
sigmah = sigma_h(x,y,curElem,lvl,p);

val = sum( (sigmah - Asigmah).*(sigmah - Asigmah),2 );
val = reshape(val,[1 1 length(x)]);
%%
function evalEnergy = getConjEnergy(x,y,curElem,lvl,p)

% conjFunctional = p.problem.conjNonLinearFunc;
conjFunctional = p.problem.nonLinearConjExact;
sigma_h = p.statics.sigma_h;

evalSigma = sigma_h(x,y,curElem,lvl,p);
absSigma = (evalSigma(:,1).^2 + evalSigma(:,2).^2).^(1/2);

evalEnergy = conjFunctional(absSigma,curElem,lvl,p);

evalEnergy = reshape(evalEnergy,[1 1 length(x)]);

% %%
% function val = getAitken(x)
% 
% val = zeros(length(x)-2,1);
% 
% for k = 1:length(x)-2
%     val(k) = x(k) - (x(k+1)-x(k))^2/...
%                     (x(k+2)-2*x(k+1)+x(k));
% end
% 
% %% 
function p = drawFigures(p,k)

if k < 3
    figure(1)
    set(gcf,'Name','Square');
    if k == 1
        p.params.output.name = 'uniform \eta_1';
        p.params.output.linestyle = '-.';
        p.params.output.marker = 's';
        p = show('drawError_estimatedError',p);
%         saveas(gcf,'Square_adaptVsUnif.fig')
    elseif k == 2
        p.params.output.name = 'adaptive \eta_1';
        p.params.output.linestyle = ':';
        p.params.output.marker = 'p';
        p = show('drawError_estimatedError',p);
%         saveas(gcf,'Square_adaptVsUnif.fig')
    end
elseif k < 5
    figure(2)
    set(gcf,'Name','Lshape');
    if k == 3
        p.params.output.name = 'uniform \eta_1';
         p.params.output.linestyle = '-.';
        p.params.output.marker = 's';
        p = show('drawError_estimatedError',p);
%         saveas(gcf,'Lshape_adaptVsUnif.fig')
    elseif k == 4
        p.params.output.name = 'adaptive \eta_1';
        p.params.output.linestyle = ':';
        p.params.output.marker = 'p';
        p = show('drawError_estimatedError',p);
%         saveas(gcf,'Lshape_adaptVsUnif.fig')
    end
 elseif k < 7
    figure(3)
    set(gcf,'Name','Octagon');
    if k == 5
        p.params.output.name = 'uniform \eta_1';
         p.params.output.linestyle = '-.';
        p.params.output.marker = 's';
        p = show('drawError_estimatedError',p);
%         saveas(gcf,'Octagon_adaptVsUnif.fig')
    elseif k == 6
        p.params.output.name = 'adaptive \eta_1';
        p = show('drawError_estimatedError',p);
        p.params.output.linestyle = ':';
        p.params.output.marker = 'p';
%         saveas(gcf,'Octagon_adaptVsUnif.fig')
    end
 elseif k < 9
    figure(4)
    set(gcf,'Name','SquareSlit');
    if k == 7
        p.params.output.name = 'uniform \eta_1';
         p.params.output.linestyle = '-.';
        p.params.output.marker = 's';
        p = show('drawError_estimatedError',p);
%         saveas(gcf,'SquareSlit_adaptVsUnif.fig')
    elseif k == 8
        p.params.output.name = 'adaptive \eta_1';
        p.params.output.linestyle = ':';
        p.params.output.marker = 'p';
        p = show('drawError_estimatedError',p);
%         saveas(gcf,'SquareSlit_adaptVsUnif.fig')
    end   
end
