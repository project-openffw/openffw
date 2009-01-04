function p = elasticity_init(p)
% initialisations for elasticity

% Copyright 2007 Jan Reininghaus, David Guenther
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


p = elasticity_setLame(p);

p.statics.H1semiError = @H1semiErrorElasticity;
p.statics.L2error = @L2errorElasticity;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function val = L2errorElasticity(pts,curElem,lvl,p)

u_h = p.statics.u_h;
u_exact = p.problem.u_exact;

exactU =  u_exact(pts,p)';

uh =  u_h(pts,curElem,lvl,p);
val(1,1,:) = sum((uh - exactU).*(uh - exactU),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = H1semiErrorElasticity(pts,curElem,lvl,p)

sigma_h = p.statics.sigma_h;
sigma_exact = p.problem.sigma_exact;

exactSigma =  sigma_exact(pts,p);
sigmah =  sigma_h(pts,curElem,lvl,p);

val(1,1,:) = sum(sum( (sigmah - exactSigma).*(sigmah - exactSigma),2 ),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
