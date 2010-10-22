function p = functionalMeanIntegral(p,curLvl)

if nargin < 2
    curLvl = length(p.level);
end

% First Goal functional is Integral mean value 
% of v over [0,1] \times [-1,1] SquareSlit
% of v over [0,1] \times [0,1] Lshape

% Second Goal functional is Integral mean value 
% of Dv over x in [-1,0.5] SquareSlit
% of Dv over x in [-1,0.5] Lshape

% Third Goal functional is L^2-Norm 
% of DW(Dv)=sigma over Omega

%% input
n4e = p.level(curLvl).geom.n4e;
c4n = p.level(curLvl).geom.c4n;
degree = 2;

%%
x = c4n(:,1);
x = x(n4e);
y = c4n(:,2);
y = y(n4e);

GOAL = p.params.GOAL;
if strcmp(GOAL,'GOAL_IM')

    ind = find(x(:,1)>= 0 & x(:,2)>= 0 & x(:,3)>= 0);

    if isempty(ind)
        value = 0;
    else
%       value = 0.5*integrate(n4e(ind,:), curLvl,  degree, @integrand, p);    %SquareSlit
        value = 1*integrate(n4e(ind,:), curLvl,  degree, @integrand, p);    %Lshape
    end

elseif strcmp(GOAL,'GOAL_DP')
    ind = find(x(:,1)<= -0.5 & x(:,2)<= -0.5 & x(:,3)<= -0.5);

    if isempty(ind)
        value = 0;
    else
        value = 1*integrate(n4e(ind,:), curLvl,  degree, @integrand, p);    %Lshape + SquareSlit
    end
else
    ind = find(x(:,1)<= 1 & x(:,2)<= 1 & x(:,3)<= 1); % all x

    if isempty(ind)
        value = 0;
    else
        value = 1*integrate(n4e(ind,:), curLvl,  degree, @integrand, p);    %Lshape + SquareSlit
    end
end

%% Output
p.level(curLvl).Ju = sum(value);


% Integral Mean Value of u_h
function val = integrand(x,y,parts,lvl,p)
u_h = p.statics.u_h;
uh =  u_h(x,y,parts,lvl,p);
val(1,1,:) = uh(:);

