function val = P3getBasisCoefficients()

% author: Joscha Gedicke

% val = [lambdaVector(0,0),lambdaVector(1,0),lambdaVector(0,1),...
%        lambdaVector(1/3,0),lambdaVector(2/3,1/3),lambdaVector(0,2/3),...
%        lambdaVector(2/3,0),lambdaVector(1/3,2/3),lambdaVector(0,1/3),...
%        lambdaVector(1/3,1/3)];
% val = val^-1;

val = [...
    1.0000         0         0   -2.2500         0   -2.2500   -2.2500         0    2.2500    4.5000
         0    1.0000         0   -2.2500   -2.2500         0    2.2500   -2.2500         0    4.5000
         0         0    1.0000         0   -2.2500   -2.2500         0    2.2500   -2.2500    4.5000
         0         0         0    2.2500         0         0    6.7500         0         0   -6.7500
         0         0         0         0    2.2500         0         0    6.7500         0   -6.7500
         0         0         0         0         0    2.2500         0         0    6.7500   -6.7500
         0         0         0    2.2500         0         0   -6.7500         0         0   -6.7500
         0         0         0         0    2.2500         0         0   -6.7500         0   -6.7500
         0         0         0         0         0    2.2500         0         0   -6.7500   -6.7500
         0         0         0         0         0         0         0         0         0   27.0000
     ];
 
function val = lambdaVector(x,y)
val = [ lambda1(x,y);
        lambda2(x,y);
        lambda3(x,y);
        lambda1(x,y)*lambda2(x,y);
        lambda2(x,y)*lambda3(x,y);
        lambda3(x,y)*lambda1(x,y);
        lambda1(x,y)*lambda2(x,y)*(lambda1(x,y)-lambda2(x,y));
        lambda2(x,y)*lambda3(x,y)*(lambda2(x,y)-lambda3(x,y));
        lambda3(x,y)*lambda1(x,y)*(lambda3(x,y)-lambda1(x,y));
        lambda1(x,y)*lambda2(x,y)*lambda3(x,y)];

function val = lambda1(x,y)
val = 1-x-y;
function val = lambda2(x,y)
val = x;
function val = lambda3(x,y)
val = y;
