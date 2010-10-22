function p = functionalMeanIntegralExact(p)

curLvl = length(p.level);

% Target functional is Integral mean value 
% of u over [0.5,1] \times [0,1]

%% input
n4e = p.level(curLvl).geom.n4e;
c4n = p.level(curLvl).geom.c4n;
degree = 1;

%%
x = c4n(:,1);
x = x(n4e);
y = c4n(:,2);
y = y(n4e);

ind = find(x(:,1)>= 0.5 & x(:,2)>= 0.5 & x(:,3)>= 0.5);
ind2 = find(x(:,1)< 0.5 | x(:,2)< 0.5 | x(:,3)< 0.5);

if isempty(ind)
    value = 0;
    value2 = 0;
else
    value = 2*integrate(n4e(ind,:), curLvl,  degree, @integrand, p);
    value2 = 0*integrate(n4e(ind2,:), curLvl,  degree, @integrand, p);
end

%% Output
p.MeanIntegralExact = sum(value) + sum(value2);


% Integral Mean Value of u_h
function val = integrand(pts,pts_ref,parts,lvl,p)
u_exact = p.problem.u_exact;
u =  u_exact(pts,pts_ref,p);
val(1,1,:) = u(:);

