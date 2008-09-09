function p = gettingStarted

problem = 'Elliptic_Lshape';
pdeSolver = 'P1';
maxNrDoF = 100;
mark = 'bulk';

p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elliptic');
p = p.statics.run(p);

figure(1);
p = show('drawU',p);

figure(2);
p = show('drawError_estimatedError',p);