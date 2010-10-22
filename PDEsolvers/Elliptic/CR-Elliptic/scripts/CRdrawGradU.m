function p = CRdrawGradU(p,lvl)
% draw piecewise constant vector field with default parameters:
% drawInfo = true
% localRes = 5
% drawStreamslice = true

if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set graphic options from structure p
drawInfo = loadField('p.params.output.drawGradU','drawInfo',p,true);
localRes = loadField('p.params.output.drawGradU','localRes',p,5);
drawStreamslice = loadField('p.params.output.drawGradU','drawStreamslice',p,true);

% load geometry
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

% load gradients
grad = p.level(lvl).grad4e;

% load enumerated data
nrElems = p.level(lvl).nrElems;
nrDoF = p.level(lvl).nrDoF;
area4e = p.level(lvl).enum.area4e;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nrElems > 50)
	warning('This may take forever. Please try a smaller level')
	pause
end

clf;
hold on;

coordX = reshape(c4n(n4e',1),3,nrElems);
coordY = reshape(c4n(n4e',2),3,nrElems);
normGrad = (grad(:,1).^2 + grad(:,2).^2).^(1/2);
patch(coordX,coordY,(normGrad*[1 1 1])');

if(drawStreamslice)
    for curElem = 1:nrElems
        curU = grad(curElem,1);
        curV = grad(curElem,2);
        curArea = area4e(curElem);
        curNodes = n4e(curElem,:);
        coordX = c4n(curNodes,1);
        coordY = c4n(curNodes,2);
        minX = min(coordX);
        minY = min(coordY);
        maxX = max(coordX);
        maxY = max(coordY);

        diam = sqrt(curArea);
        XI = minX:diam/localRes:maxX;
        YI = minY:diam/localRes:maxY;

        UI=tri2grid([coordX,coordY]',[1,2,3,1]' ,...
            [curU;curU;curU],XI,YI' );
        VI=tri2grid([coordX,coordY]',[1,2,3,1]' ,...
            [curV;curV;curV],XI,YI' );

        % Draw Streamlines
        handle = streamslice(XI,YI,UI,VI);
        set(handle,'Color','black');
    end
end

hold off

if(drawInfo)
	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
	title('Piecewise constant gradient of CR Solution');
end

