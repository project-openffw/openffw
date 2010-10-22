function plotDummy
clear all
index = {'1','01','001'};

for k = 1
    string = ['optimalBetaFactor' '1' 'SquareExact.mat'];
    load(string)
    determOptimalParam2(A,'1')
end


%%
function A = determOptimalParam2(A,val)

nrStructs = length(A.struct);

for k = 1:nrStructs
    p = A.struct(k).p;
    p.params.output.fontSize = 18;
    p.params.output.holdIt = true;
%     for lvl = 2:length(p.level)
%         [k,lvl]
%         conjEnergyExact = integrate(p.level(lvl).geom.n4e,lvl,10,@getConjEnergy,p);
%         p.level(lvl).conjEnergyExact = conjEnergyExact;
%     end
%     for lvl = 2:length(p.level)
% %         [k,lvl]
%         dof(lvl-1) = p.level(lvl).nrDoF;
%         energy(lvl-1) = sum(p.level(lvl).conjEnergyExact);
%     end
%     p.level(end).refineEdges = true(p.level(end).nrEdges,1);
%     p = redGreenBlue(p);
%     p = ODRTenumerate(p);
%     level = length(p.level);
%     limit = sum(integrate(p.level(level).geom.n4e,level,10,@getConjEnergyExact,p));
%     figure(1)
%     p = show('drawError_estimatedError',p);
%     string1 = ['optimalBetaFactor' val 'SquareExactEta1'];
%     saveas(gcf,string1,'epsc2');
%     saveas(gcf,string1,'fig');
    
%     figure(2)
    q = p;
    for k = 2:length(p.level)
        w.problem = p.problem;
        w.params = p.params;
        w.statics = p.statics;
        w.level = p.level(k);
        w = ODRTestimate_Proj(w);
        p.level(k).etaProj = w.level(end).estimatedError;
        q.level(k).estimatedError = p.level(k).etaProj;
        q.level(k).etaProj4e = w.level(end).etaT;
    end
%     q = show('drawError_estimatedError',q);
%     string2 = ['optimalBetaFactor' val 'SquareExactEta2'];
%     saveas(gcf,string2,'epsc2');
%     saveas(gcf,string2,'fig');
%     
%     figure(3)
%     p = show('drawError_L2errorPhminusP0',p);
%     string3 = ['optimalBetaFactor' val 'SquareExactStress'];
%     saveas(gcf,string3,'epsc2');
%     saveas(gcf,string3,'fig');
    
%     A.struct(k).p = p;
%     origEnergy = energy
%     aitken = getAitken(energy)
%     limit = aitken(end)
%     errorEnergy = abs(origEnergy-limit)
%     loglog(dof,errorEnergy);
%     hold all
%     minDoF = 1000;
%     I = find(dof >= minDoF);
%     nrDoF4lvl = dof(I);
%     err4lvl = errorEnergy(I);
% 
%     nrDoF4lvl = log(nrDoF4lvl);
%     err4lvl = log(err4lvl);
%     approx = polyfit(nrDoF4lvl,err4lvl,1);
% 
%     y = polyval(approx,nrDoF4lvl);
% %     rate = -(y(end) - y(1)) / (nrDoF4lvl(end) - nrDoF4lvl(1))
    volumeFraction = p.statics.volumeFraction;
%     q = p;
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
%         I = find((dummy == 3) | (dummy == 0));
%         I = find((dummy == 3));
        I = find((dummy == 0));
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

        etaProj = q.level(lvl).etaProj4e;
        innerError = etaProj(I);
        
        estimatedError = sqrt(sum(innerError.^2));
        q.level(lvl).estimatedError = estimatedError;
    end
    q = show('drawError_estimatedError',q);
    clear p q;
end

% string = ['optimalBetaFactor' val 'SquareExact.mat'];
% save(string,'A');

%%
function evalEnergy = getConjEnergy(x,y,curElem,lvl,p)

% conjFunctional = p.problem.conjNonLinearFunc;
conjFunctional = p.problem.nonLinearConjExact;
sigma_h = p.statics.sigma_h;

evalSigma = sigma_h(x,y,curElem,lvl,p);
absSigma = (evalSigma(:,1).^2 + evalSigma(:,2).^2).^(1/2);

evalEnergy = conjFunctional(absSigma,curElem,lvl,p);

evalEnergy = reshape(evalEnergy,[1 1 length(x)]);

%%
function evalEnergy = getConjEnergyExact(x,y,curElem,lvl,p)

% conjFunctional = p.problem.conjNonLinearFunc;
conjFunctional = p.problem.nonLinearConjExact;
sigma = p.problem.sigma0;

evalSigma = sigma(x,y,curElem,lvl,p);
absSigma = (evalSigma(:,1).^2 + evalSigma(:,2).^2).^(1/2);

evalEnergy = conjFunctional(absSigma,curElem,lvl,p);

evalEnergy = reshape(evalEnergy,[1 1 length(x)]);

%%
function val = getAitken(x)

val = zeros(length(x)-2,1);

for k = 1:length(x)-2
    val(k) = x(k) - (x(k+1)-x(k))^2/...
                    (x(k+2)-2*x(k+1)+x(k));
end