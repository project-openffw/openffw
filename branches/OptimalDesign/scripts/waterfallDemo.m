function p = waterfallDemo

warning off
restoredefaultpath
warning on
addpath(fullfile(pwd,'scripts'));
addpath(fullfile(pwd,'helpers'));

problem = 'Laplace-Square-exact';

% pdeSolver = 'CR-Elliptic';
% pdeSolver = 'RT0-P0-mixed-Elliptic';
pdeSolver = 'P1-Elliptic';

maxNrDoF = 3900;

mark = 'bulk'
% mark = 'uniform'

configuration = [];

k = [1:10:51];
level = 1;
clf
for curK = k
	curK
	clear p
	%%%%%%%%%%%%%%%%%%%%%%COMPUTE DISCRETE SOLUTION%%%%%%%%%%%%%%%%%%%%%%%
	mark = 'bulk'
	p = configureP(pdeSolver,problem,mark,maxNrDoF,configuration);
	p.problem.k = curK;
	p = run(p);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	
	% saveFigures = true;
	saveFigures = false;
	%format = '.fig';
	format = '.jpg';
	%format = '.eps';

	% drawFunction = 'drawU';
	% drawFunction = 'drawGrid';

	% drawFunction = 'drawError_approxL2error';
	% drawFunction = 'drawError_estimatedError';
	drawFunction = 'drawError_approxH1semiError';

	% clf
% 	hold all
	p = show(drawFunction,saveFigures,format,p);
	% hold off

	nrDoF = p.level(end).nrDoF;
	bulkError = p.level(end).approxH1semiError;
	bulkError = bulkError/sqrt(nrDoF);


	% pause;
	% drawFunction = 'drawU';
	% p = show(drawFunction,saveFigures,format,p);


	mark = 'uniform'
	%%%%%%%%%%%%%%%%%%%%%%COMPUTE DISCRETE SOLUTION%%%%%%%%%%%%%%%%%%%%%%%
	% drawFunction = 'drawError_approxH1semiError';
	p = configureP(pdeSolver,problem,mark,maxNrDoF,configuration);
	p.problem.k = curK;
	p = run(p);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	hold all
	p = show(drawFunction,saveFigures,format,p);
% 	hold off

	nrDoF = p.level(end).nrDoF;
	unifError = p.level(end).approxH1semiError;
	unifError = unifError/sqrt(nrDoF);

	quotient(level) = unifError/bulkError
	level = level +1;
end

figure(2)
plot(k,quotient)

function p = configureP(pdeSolver,problem,mark,maxNrDoF,configuration)
disp(pdeSolver)
p = loadConfig(configuration);
p.params.mark = mark;
p.params.problem.name = problem;
p.params.pdeSolver = pdeSolver;
p.params.maxNrDoF = maxNrDoF;
p.params.problem.type = 'elliptic';
p = initFFW(p);
