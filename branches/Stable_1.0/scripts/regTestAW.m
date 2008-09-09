function p = startElasticity

warning off
restoredefaultpath
warning on
addpath(fullfile(pwd,'scripts'));
addpath(fullfile(pwd,'helpers'));

problem = 'Elasticity-Lshape';
% problem = 'Elasticity-Cooks';

maxNrDoF = 1000;
% mark = 'bulk'
mark = 'uniform'


pdeSolver = 'Arnold-Winther-mixed-Elasticity-Paper';


%%%%%%%%%%%%%%%%%%%%%%COMPUTE DISCRETE SOLUTION%%%%%%%%%%%%%%%%%%%%%%%
p = configureP(pdeSolver,problem,mark,maxNrDoF,[]);
p = run(p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xNew = p.level(end).x;


pdeSolver = 'Arnold-Winther-mixed-Elasticity';
% %%%%%%%%%%%%%%%%%%%%%%COMPUTE DISCRETE SOLUTION%%%%%%%%%%%%%%%%%%%%%%%
p = configureP(pdeSolver,problem,mark,maxNrDoF,[]);
p = run(p);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xOld = p.level(end).x;

Error = norm(xOld-xNew)



function p = configureP(pdeSolver,problem,mark,maxNrDoF,configuration)
	disp(pdeSolver)
	p = loadConfig(configuration);
	p.params.mark = mark;
	p.params.problem.name = problem;
	p.params.pdeSolver = pdeSolver;
	p.params.maxNrDoF = maxNrDoF;
	p.params.problem.type = 'elasticity';
	p = initFFW(p);










