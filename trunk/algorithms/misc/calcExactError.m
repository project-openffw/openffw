function p = calcExactError(errType,curLvl,p)
% calcExactError return the exact error to the discrete solution in a norm
% dependent on 'errType', i.e., H1semi or L2.

% Copyright 2007 Andreas Byfut, Joscha Gedicke
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

%% calculate exact error
degree = p.params.integrationDegrees.exactError;
intMode = loadField('p.params','IntegrationMode',p,'elementwise');

if isempty(loadField('p.problem','u_exact',p))
    error('Cannot compute exact error because no exact solution was specified');
end

fprintf('\t Integrating, Level %2i of %2i, nrGaussPoints= %i',curLvl,p.level(end).level...
    , ceil((degree+1)/2 )^2 );
time = cputime;
if strcmpi(intMode,'elementwise')
    eval(['integrandFunc = p.statics.',errType,';']);
    err4e = integrate( p.level(curLvl).geom.n4e, curLvl,  degree, integrandFunc, p );
else
    eval(['integrandFunc = p.statics.',errType,'Vectorised;']);
    err4e = integrateVectorised( p.level(curLvl).geom.n4e, curLvl,  degree, integrandFunc, p );
end
time = cputime - time;
fprintf(', elapsed time = %.1g sec \n',time);
eval(['p.level(curLvl).',errType,' = sqrt(sum(err4e));']);
eval(['p.level(curLvl).',errType,'4e = sqrt(err4e);']);
return