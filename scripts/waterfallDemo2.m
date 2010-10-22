function p = waterfallDemo(p)

warning off
restoredefaultpath
warning on
addpath(fullfile(pwd,'scripts'));
addpath(fullfile(pwd,'helpers'));
addpath(fullfile(pwd,'algorithms','misc'));

problem = 'Laplace-Square-exact';


pdeSolver = 'P1-Elliptic';
% pdeSolver = 'CR-Elliptic';
% pdeSolver = 'RT0-P0-mixed-Elliptic';

maxNrDoF = 3000;

mark = 'bulk'
% mark = 'uniform'

%%%%%%%%%%%%%%%%%%%%%COMPUTE DISCRETE SOLUTION%%%%%%%%%%%%%%%%%%%%%%%
p = configureP(pdeSolver,problem,mark,maxNrDoF,[],'elliptic');
p.problem.k = 60;
p = run(p);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); hold all; p = show('drawError_H1semiError',p);
figure(2); hold off; p = show('drawU',p);

%%%%%%%%%%%%%%%%%%%%%COMPUTE DISCRETE SOLUTION%%%%%%%%%%%%%%%%%%%%%%%
p = configureP(pdeSolver,problem,mark,maxNrDoF,[],'elliptic');
p.problem.k = 15;
p = run(p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); hold all; p = show('drawError_H1semiError',p);
figure(3); hold off; p = show('drawU',p);






