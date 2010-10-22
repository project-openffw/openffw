function drawDisplacementBasis(j,res,p)

u = p.statics.basisU;

XI = 0:1/res:1;
YI = 0:1/res:1;

UI = zeros(res);

for indX = 1:length(XI)
	for indY = 1:length(YI)
        x = XI(indX);
        y = YI(indY);
        if(x+y > 1)
            val = [NaN,NaN,NaN];
        else
            val  = u(x,y,1,1,p);
        end
        UI(indX,indY) = val(j);
	end
end

surf(XI,YI,UI','EdgeColor','none');