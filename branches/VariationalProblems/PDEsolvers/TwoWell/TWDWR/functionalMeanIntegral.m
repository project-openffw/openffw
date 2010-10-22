function p = functionalMeanIntegral(p,curLvl)

if nargin < 2
    curLvl = length(p.level);
end

% Target functional is Integral mean value 
% of u over [0,1] \times [-1,1] SquareLit
% of u over [0,1] \times [0,1] Lshape

%% input
n4e = p.level(curLvl).geom.n4e;
c4n = p.level(curLvl).geom.c4n;
degree = 2;

%%
x = c4n(:,1);
x = x(n4e);
y = c4n(:,2);
y = y(n4e);

ind = find(x(:,1)>= 0 & x(:,2)>= 0 & x(:,3)>= 0);

if isempty(ind)
    value = 0;
else
%    value = 0.5*integrate(n4e(ind,:), curLvl,  degree, @integrand, p);    %SquareSlit
    value = 1*integrate(n4e(ind,:), curLvl,  degree, @integrand, p);    %Lshape
end

%% Output
p.level(curLvl).Ju = sum(value);


% Integral Mean Value of u_h
function val = integrand(x,y,parts,lvl,p)
u_h = p.statics.u_h;
uh =  u_h(x,y,parts,lvl,p);
val(1,1,:) = uh(:);

