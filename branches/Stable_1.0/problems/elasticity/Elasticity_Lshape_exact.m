function p = Elasticity_Lshape_exact(p)

% PDE definition
p.problem.geom = 'Lshape3';
p.problem.f = @f;
p.problem.g = @g;
p.problem.u_D = @u_D;
p.problem.u_exact = @u_exact;
p.problem.alph = .544483737;
p.PDE.E = 100000; 
p.PDE.nu = 0.3;


p.problem.sigma_exact = @sigma;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Volume force
function z = f(x,y,p)
z = zeros(size(x,1),2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dirichlet boundary values
function z = u_D(x,y,p)

z = p.problem.u_exact(x,y,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exact sigma = C*eps(u)
function z = sigma(x,y,p)

warning off;
if (x(1) >= 0 && y(1) >= 0)
    z = ([[-(cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*x.^2-cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x)).*p.problem.alph.*x.^2+3.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x)).*x.^2+cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*x.^2+2.*x.*y.*cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x))+2.*x.*y.*sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph+2.*x.*y.*sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi)-2.*x.*y.*cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*p.problem.alph+cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x)).*y.^2-cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*y.^2+cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x)).*p.problem.alph.*y.^2-cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*y.^2).*(x.^2+y.^2).^(1./2.*p.problem.alph-3./2).*p.problem.alph.*sign(x)./cos(3./4.*(p.problem.alph-1).*pi)+x-x+y-y,-(-cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*x.^2+cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*p.problem.alph.*x.^2-sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*x.^2-sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*x.^2-2.*x.*y.*p.problem.alph.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x))+2.*x.*y.*p.problem.alph.*cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi)+2.*x.*y.*cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi)+2.*x.*y.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x))+cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*y.^2-cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*p.problem.alph.*y.^2+sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*y.^2+sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*y.^2).*p.problem.alph.*(x.^2+y.^2).^(1./2.*p.problem.alph-3./2).*sign(x)./cos(3./4.*(p.problem.alph-1).*pi)+x-x+y-y];[-(-cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*x.^2+cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*p.problem.alph.*x.^2-sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*x.^2-sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*x.^2-2.*x.*y.*p.problem.alph.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x))+2.*x.*y.*p.problem.alph.*cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi)+2.*x.*y.*cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi)+2.*x.*y.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x))+cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*y.^2-cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*p.problem.alph.*y.^2+sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*y.^2+sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*y.^2).*p.problem.alph.*(x.^2+y.^2).^(1./2.*p.problem.alph-3./2).*sign(x)./cos(3./4.*(p.problem.alph-1).*pi)+x-x+y-y,(x.^2+y.^2).^(1./2.*p.problem.alph-3./2).*p.problem.alph.*(-cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x)).*x.^2+cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*x.^2-cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x)).*p.problem.alph.*x.^2+cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*x.^2+2.*x.*y.*sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph+2.*x.*y.*cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x))-2.*x.*y.*cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*p.problem.alph+2.*x.*y.*sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi)-3.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x)).*y.^2-cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*y.^2+cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x)).*p.problem.alph.*y.^2-cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*y.^2)./cos(3./4.*(p.problem.alph-1).*pi).*sign(x)+x-x+y-y]]);
elseif (x(1) <= 0 && y(1) >= 0)
	z = ([[(x.^2+y.^2).^(1./2.*p.problem.alph-3./2).*p.problem.alph.*(cos((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*x.^2+cos((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi).*x.^2-cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(pi+atan(y./x))).*p.problem.alph.*x.^2+3.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(pi+atan(y./x))).*x.^2-2.*x.*y.*cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(pi+atan(y./x))).*p.problem.alph+2.*x.*y.*sin((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi)+2.*x.*y.*sin((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph+2.*x.*y.*cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(pi+atan(y./x)))+cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(pi+atan(y./x))).*y.^2-cos((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*y.^2-cos((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi).*y.^2+cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(pi+atan(y./x))).*p.problem.alph.*y.^2)./cos(3./4.*(p.problem.alph-1).*pi).*sign(x)+x-x+y-y,-(-cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(pi+atan(y./x))).*p.problem.alph.*x.^2+sin((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*x.^2+sin((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi).*x.^2+cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(pi+atan(y./x))).*x.^2-2.*x.*y.*p.problem.alph.*cos((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi)-2.*x.*y.*cos((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi)-2.*x.*y.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(pi+atan(y./x)))+2.*x.*y.*p.problem.alph.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(pi+atan(y./x)))+cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(pi+atan(y./x))).*p.problem.alph.*y.^2-sin((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*y.^2-sin((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi).*y.^2-cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(pi+atan(y./x))).*y.^2).*p.problem.alph.*(x.^2+y.^2).^(1./2.*p.problem.alph-3./2).*sign(x)./cos(3./4.*(p.problem.alph-1).*pi)+x-x+y-y];[-(-cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(pi+atan(y./x))).*p.problem.alph.*x.^2+sin((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*x.^2+sin((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi).*x.^2+cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(pi+atan(y./x))).*x.^2-2.*x.*y.*p.problem.alph.*cos((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi)-2.*x.*y.*cos((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi)-2.*x.*y.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(pi+atan(y./x)))+2.*x.*y.*p.problem.alph.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(pi+atan(y./x)))+cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(pi+atan(y./x))).*p.problem.alph.*y.^2-sin((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*y.^2-sin((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi).*y.^2-cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(pi+atan(y./x))).*y.^2).*p.problem.alph.*(x.^2+y.^2).^(1./2.*p.problem.alph-3./2).*sign(x)./cos(3./4.*(p.problem.alph-1).*pi)+x-x+y-y,(x.^2+y.^2).^(1./2.*p.problem.alph-3./2).*p.problem.alph.*(cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(pi+atan(y./x))).*p.problem.alph.*x.^2-cos((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*x.^2-cos((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi).*x.^2+cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(pi+atan(y./x))).*x.^2+2.*x.*y.*cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(pi+atan(y./x))).*p.problem.alph-2.*x.*y.*cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(pi+atan(y./x)))-2.*x.*y.*sin((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi)-2.*x.*y.*sin((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph+3.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(pi+atan(y./x))).*y.^2-cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(pi+atan(y./x))).*p.problem.alph.*y.^2+cos((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*y.^2+cos((p.problem.alph+1).*(pi+atan(y./x))).*cos(3./4.*(p.problem.alph-1).*pi).*y.^2)./cos(3./4.*(p.problem.alph-1).*pi).*sign(x)+x-x+y-y]]);
elseif (x(1) <= 0 && y(1) <= 0)
	z = ([[-(x.^2+y.^2).^(1./2.*p.problem.alph-3./2).*p.problem.alph.*(-cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(atan(y./x)-pi)).*y.^2+cos((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*y.^2-cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(atan(y./x)-pi)).*p.problem.alph.*y.^2+cos((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*y.^2-2.*y.*x.*sin((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph+2.*y.*x.*cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(atan(y./x)-pi)).*p.problem.alph-2.*y.*x.*cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(atan(y./x)-pi))-2.*y.*x.*sin((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi)-3.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(atan(y./x)-pi)).*x.^2-cos((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*x.^2+cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(atan(y./x)-pi)).*p.problem.alph.*x.^2-cos((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*x.^2)./cos(3./4.*(p.problem.alph-1).*pi).*sign(x)+x-x+y-y,-(cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(atan(y./x)-pi)).*p.problem.alph.*y.^2-sin((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*y.^2-cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(atan(y./x)-pi)).*y.^2-sin((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*y.^2-2.*y.*p.problem.alph.*x.*cos((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi)-2.*y.*cos((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*x-2.*y.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(atan(y./x)-pi)).*x+2.*y.*p.problem.alph.*x.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(atan(y./x)-pi))-cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(atan(y./x)-pi)).*p.problem.alph.*x.^2+sin((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*x.^2+cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(atan(y./x)-pi)).*x.^2+sin((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*x.^2).*p.problem.alph.*(x.^2+y.^2).^(1./2.*p.problem.alph-3./2).*sign(x)./cos(3./4.*(p.problem.alph-1).*pi)+x-x+y-y];[-(cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(atan(y./x)-pi)).*p.problem.alph.*y.^2-sin((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*y.^2-cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(atan(y./x)-pi)).*y.^2-sin((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*y.^2-2.*y.*p.problem.alph.*x.*cos((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi)-2.*y.*cos((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*x-2.*y.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(atan(y./x)-pi)).*x+2.*y.*p.problem.alph.*x.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(atan(y./x)-pi))-cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(atan(y./x)-pi)).*p.problem.alph.*x.^2+sin((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*x.^2+cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(atan(y./x)-pi)).*x.^2+sin((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*x.^2).*p.problem.alph.*(x.^2+y.^2).^(1./2.*p.problem.alph-3./2).*sign(x)./cos(3./4.*(p.problem.alph-1).*pi)+x-x+y-y,(x.^2+y.^2).^(1./2.*p.problem.alph-3./2).*p.problem.alph.*(3.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(atan(y./x)-pi)).*y.^2+cos((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*y.^2+cos((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*y.^2-cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(atan(y./x)-pi)).*p.problem.alph.*y.^2-2.*y.*x.*sin((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi)+2.*y.*x.*cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(atan(y./x)-pi)).*p.problem.alph-2.*y.*x.*cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*(atan(y./x)-pi))-2.*y.*x.*sin((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph+cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(atan(y./x)-pi)).*x.^2-cos((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*x.^2-cos((p.problem.alph+1).*(atan(y./x)-pi)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*x.^2+cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*(atan(y./x)-pi)).*p.problem.alph.*x.^2)./cos(3./4.*(p.problem.alph-1).*pi).*sign(x)+x-x+y-y]]);
elseif (x(1) >= 0 && y(1) <= 0)
	z = ([[-(x.^2+y.^2).^(1./2.*p.problem.alph-3./2).*p.problem.alph.*(3.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x)).*x.^2-cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x)).*p.problem.alph.*x.^2+cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*x.^2+cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*x.^2+2.*x.*y.*sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph-2.*x.*y.*cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*p.problem.alph+2.*x.*y.*cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x))+2.*x.*y.*sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi)+cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x)).*y.^2+cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x)).*p.problem.alph.*y.^2-cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*y.^2-cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*y.^2).*sign(x)./cos(3./4.*(p.problem.alph-1).*pi)+x-x+y-y,-(-sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*x.^2-cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*x.^2+cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*p.problem.alph.*x.^2-sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*x.^2+2.*x.*y.*p.problem.alph.*cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi)+2.*x.*y.*cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi)-2.*x.*y.*p.problem.alph.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x))+2.*x.*y.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x))+sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*y.^2+cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*y.^2-cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*p.problem.alph.*y.^2+sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*y.^2).*p.problem.alph.*(x.^2+y.^2).^(1./2.*p.problem.alph-3./2).*sign(x)./cos(3./4.*(p.problem.alph-1).*pi)+x-x+y-y];[-(-sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*x.^2-cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*x.^2+cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*p.problem.alph.*x.^2-sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*x.^2+2.*x.*y.*p.problem.alph.*cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi)+2.*x.*y.*cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi)-2.*x.*y.*p.problem.alph.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x))+2.*x.*y.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x))+sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*y.^2+cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*y.^2-cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*p.problem.alph.*y.^2+sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*y.^2).*p.problem.alph.*(x.^2+y.^2).^(1./2.*p.problem.alph-3./2).*sign(x)./cos(3./4.*(p.problem.alph-1).*pi)+x-x+y-y,-(x.^2+y.^2).^(1./2.*p.problem.alph-3./2).*p.problem.alph.*(cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x)).*x.^2+cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x)).*p.problem.alph.*x.^2-cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*x.^2-cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*x.^2-2.*x.*y.*sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph+2.*x.*y.*cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x)).*p.problem.alph-2.*x.*y.*cos(3./4.*(p.problem.alph+1).*pi).*sin((p.problem.alph-1).*atan(y./x))-2.*x.*y.*sin((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi)+3.*cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x)).*y.^2-cos(3./4.*(p.problem.alph+1).*pi).*cos((p.problem.alph-1).*atan(y./x)).*p.problem.alph.*y.^2+cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*y.^2+cos((p.problem.alph+1).*atan(y./x)).*cos(3./4.*(p.problem.alph-1).*pi).*p.problem.alph.*y.^2).*sign(x)./cos(3./4.*(p.problem.alph-1).*pi)+x-x+y-y]]);
end
warning on;

if(ndims(z) < 3)
	z = reshape(z,[],4)';
	z = reshape(z,2,2,[]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%exact solution
function z = u_exact(x,y,p)

mu = p.PDE.mu;
lambda = p.PDE.lambda;

[phi,r] = cart2pol(x,y);
alph = p.problem.alph;
omega = 3*pi/4;

C_1 = -cos((alph+1)*omega)/cos((alph-1)*omega);
C_2 = 2*(lambda+2*mu)/(lambda+mu);

ut = (1/(2*mu)) * r.^alph .*((alph+1)*sin((alph+1)*phi)+ ...
    (C_2+alph-1)*C_1*sin((alph-1)*phi));
ur = (1/(2*mu))*r.^alph .* (-(alph+1)*cos((alph+1)*phi)+ ...
    (C_2-(alph+1))*C_1*cos((alph-1)*phi));
z = [ur .* cos(phi) - ut .* sin(phi), ur .* sin(phi) + ut .* cos(phi)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann boundary values
function z = g(x,y,n,p)
sigma = p.problem.sigma_exact(x,y,p);
n = reshape(n',2,1,size(n,1));
z = matMul(sigma,n);
z = squeeze(z)';
