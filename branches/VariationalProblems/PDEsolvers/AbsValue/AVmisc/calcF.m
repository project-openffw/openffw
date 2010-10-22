clear all
syms r theta r_osc eps_osc eps2_osc k real


u = r.^(2/3).*sin(2/3*theta);
% u = r.^(2/3).*sin(2/3*theta) + eps2_osc*cos(theta)*sin(eps_osc+k*r_osc)...
%         .*((r_osc/eps_osc).^4-2*(r_osc/eps_osc).^2+1);

% Volume force
matrix = [cos(theta) -1/r*sin(theta);
          sin(theta) 1/r*cos(theta)];
      
gradU = matrix*[diff(u,r); diff(u,theta)]
KappaGradU = gradU;
DiffKappaGradU = [diff(KappaGradU(1),r) diff(KappaGradU(2),r);
                  diff(KappaGradU(1),theta) diff(KappaGradU(2),theta)]
DiffKappaGradU = matrix*DiffKappaGradU;  
minusDivKappaGradU = -simple(DiffKappaGradU(1,1) + DiffKappaGradU(2,2));



f = minusDivKappaGradU 
