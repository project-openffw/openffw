function val = P3getBasisCoefficients()
% Matrix containing basis coefficients

% Copyright 2007 Joscha Gedicke
%
% This file is part of FFW.
%
% FFW is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% FFW is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


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
