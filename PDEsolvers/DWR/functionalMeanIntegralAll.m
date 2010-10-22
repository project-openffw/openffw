function p = functionalMeanIntegral(p,curLvl)

if nargin < 2
    curLvl = length(p.level);
end

%optimalDesign/AbsValue
  % First Goal functional is Integral mean value 
    % of v over [0,1]x[-1,1] Square(Slit)           1/|O1| = 0.5
    % of v over [0,1]x[0,1] Lshape                  1/|O1| = 1
  % Second Goal functional is Integral mean value 
    % of Dv over x in [-1,-0.5] Square(Slit)        1/|O1| = 1
    % of Dv over x in [-1,-0.5] Lshape              1/|O1| = 1
  % Third Goal functional is L^2-Norm 
    % of DW(Dv)=sigma over Omega                    1/|O1| = 0.25  Square(Slit) 
    %                                               1/|O1| = 1/3   Lshape

%TwoWell
  % First Goal functional is Integral mean value 
    % of v over [0.5,1]x[0,1.5] Square              1/|O1| = 4/3
  % Second Goal functional is Integral mean value 
    % of Dv over x in [0,0.25] Square               1/|O1| = 8/3
  % Third Goal functional is L^2-Norm 
    % of DW(Dv)=sigma over Omega=[0,1]x[0,1.5]      1/|O1| = 2/3

%% input
n4e = p.level(curLvl).geom.n4e;
c4n = p.level(curLvl).geom.c4n;
degree = 2;

%%
x = c4n(:,1);
x = x(n4e);
y = c4n(:,2);
y = y(n4e);

GOAL = p.params.GOAL;               %GOAL_IM, GOAL_DP, GOAL_PR
probClass = p.problem.probClass;    %optimalDesign, TwoWell, AbsValue
problem = p.params.problem.name;    %pC _ Square, SquareSlit, Lshape (_exact)
if strcmp(GOAL,'GOAL_IM')

    if (strcmp(probClass,'optimalDesign') || strcmp(probClass,'AbsValue'))
      ind = find(x(:,1)>= 0.5 & x(:,2)>= 0.5 & x(:,3)>= 0.5);
    else
      ind = find(x(:,1)>= 0.5 & x(:,2)>= 0.5 & x(:,3)>= 0.5);
    end
    
    if isempty(ind)
        value = 0;
    else
        
        if (strcmp(problem,'OptimalDesign_SquareSlit') || strcmp(problem,'OptimalDesign_SquareSlit_exact') ||...
            strcmp(problem,'AbsValue_Square') || strcmp(problem,'AbsValue_Square_exact') )
            AreaPartOmega = 0.5;
        elseif (strcmp(problem,'OptimalDesign_Lshape') || strcmp(problem,'OptimalDesign_Lshape_exact') ||...
            strcmp(problem,'AbsValue_Lshape') || strcmp(problem,'AbsValue_Lshape_exact') )
            AreaPartOmega = 1;
        else %TwoWell_Square or TwoWell_Square_exact
            AreaPartOmega = 4/3;
        end

        value = AreaPartOmega*integrate(n4e(ind,:), curLvl,  degree, @integrandIM, p);  
    end

elseif strcmp(GOAL,'GOAL_DP')

    if (strcmp(probClass,'optimalDesign') || strcmp(probClass,'AbsValue'))
      ind = find(x(:,1)>= -0.5 & x(:,2)>= -0.5 & x(:,3)>= -0.5);
    else
      ind = find(x(:,1)>= 0.25 & x(:,2)>= 0.25 & x(:,3)>= 0.25);
    end

    if isempty(ind)
        value = 0;
    else
        if (strcmp(problem,'TwoWell_Square') || strcmp(problem,'TwoWell_Square_exact') )
            AreaPartOmega = 8/3;
        else %OD + AV
            AreaPartOmega = 1;
        end

        value = AreaPartOmega*integrate(n4e(ind,:), curLvl,  degree, @integrandDP, p);  
    end
else %GOAL_PR
    ind = find(x(:,1)<= 1 & x(:,2)<= 1 & x(:,3)<= 1); % all x

    if isempty(ind)
        value = 0;
    else
        if (strcmp(problem,'OptimalDesign_SquareSlit') || strcmp(problem,'OptimalDesign_SquareSlit_exact') ||...
            strcmp(problem,'AbsValue_Square') || strcmp(problem,'AbsValue_Square_exact') )
            AreaPartOmega = 1;%0.25;
        elseif (strcmp(problem,'OptimalDesign_Lshape') || strcmp(problem,'OptimalDesign_Lshape_exact') ||...
            strcmp(problem,'AbsValue_Lshape') || strcmp(problem,'AbsValue_Lshape_exact') )
            AreaPartOmega = 1;%1/3;
        else %TwoWell_Square or TwoWell_Square_exact
            AreaPartOmega = 1;%2/3;
        end

%        if strcmp(probClass,'TwoWell')
%            value = AreaPartOmega*integrate(n4e(ind,:), curLvl,  degree, @integrandPR, p).^{3/4};  
%        else
            intTerm = integrate(n4e(ind,:), curLvl,  degree, @integrandPR, p);
            value = AreaPartOmega*intTerm;  
%        end
    end
end

%% Output
fprintf('\nJu = %.15g',sum(value))

p.level(curLvl).Ju = sum(value);


% Integral Mean Value of u_h
function val = integrandIM(x,y,parts,lvl,p)
u_h = p.statics.u_h;
uh =  u_h(x,y,parts,lvl,p);
val(1,1,:) = uh(:);

% Integral Mean Value of Du_h
function val = integrandDP(x,y,parts,lvl,p)
Du_h = p.statics.grad_h;
Duh =  Du_h(x,y,parts,lvl,p);
val(1,1,:) = sum(Duh(:),1);

% L^2-Norm of sigma_h
function val = integrandPR(x,y,parts,lvl,p)
%probClass = p.problem.probClass; 
sigma_h = p.statics.sigma_h;
evalSigma = sigma_h(x,y,parts,lvl,p);

%if strcmp(probClass,'TwoWell')
%val(1,1,:) = (evalSigma(:,1).^2 + evalSigma(:,2).^2).^{2/3}; %(x,y)^(4/3) = (x^2+y^2)^(2/3)
%else
val(1,1,:) = evalSigma(:,1).^2 + evalSigma(:,2).^2;
%end

