clear all
syms x y real


u = sin(x^3)*cos(y^3);
% u = r.^(2/3).*sin(2/3*theta) + eps2_osc*cos(theta)*sin(eps_osc+k*r_osc)...
%         .*((r_osc/eps_osc).^4-2*(r_osc/eps_osc).^2+1);

% Volume force
matrix = [2 0; 0 2];
      
gradU = matrix*[diff(u,x); diff(u,y)]
KappaGradU = gradU;
DiffKappaGradU = [diff(KappaGradU(1),x) diff(KappaGradU(2),x);
                  diff(KappaGradU(1),y) diff(KappaGradU(2),y)]
DiffKappaGradU = matrix*DiffKappaGradU;  
minusDivKappaGradU = -simple(DiffKappaGradU(1,1) + DiffKappaGradU(2,2));



f = minusDivKappaGradU;
f = Matlab4Maple(f)
