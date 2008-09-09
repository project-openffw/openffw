syms x y k 

r = sqrt((x-5/4)^2 + (y+1/4)^2);
z = x.*y.*(1-x).*(1-y).*atan(k*(r-1));
laplaceZ = simple(diff(diff(z,x),x) + diff(diff(z,y),y));
laplaceZ = char(laplaceZ);

laplaceZ = strrep(laplaceZ,'*','.*');
laplaceZ = strrep(laplaceZ,'/','./');
laplaceZ = strrep(laplaceZ,'^','.^');
laplaceZ = strrep(laplaceZ,'k','p.problem.k');

exec = ['p.problem.f = @(x,y,p)(',laplaceZ,');']
eval(exec,'disp(''err'')');
