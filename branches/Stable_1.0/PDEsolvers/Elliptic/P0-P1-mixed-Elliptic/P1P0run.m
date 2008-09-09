function p = P1P0run(p)

% Works only with u_D = 0 and no Neumann data
% (i.e. it works only on 'Sqaure' and 'Lshape')


	warning('Works only with u_D = 0 and no Neumann data');



% Initial enumeration
enumerate = p.statics.enumerate;
p = enumerate(p);

% Set initial values
nrNodes = p.level(1).nrNodes;
nrElems = p.level(1).nrElems;
p.level(1).x = zeros(2*nrNodes + nrElems,1);

p = computeSolution(p);
