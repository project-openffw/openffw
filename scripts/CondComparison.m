
warning off
restoredefaultpath
warning on
addpath(fullfile(pwd,'scripts'));
addpath(fullfile(pwd,'helpers'));



hold all;


% % %%%%%%%%%%%%%%%%%%%%%%COMPUTE DISCRETE SOLUTION%%%%%%%%%%%%%%%%%%%%%%%
p = configureP(pdeSolver,problem,mark,maxNrDoF,[]);
p.PDE.nu = nu1;
p.params.showIteratively = false;
p = run(p);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = show('drawCondNr_-',saveFigures,format,p);


% % %%%%%%%%%%%%%%%%%%%%%%COMPUTE DISCRETE SOLUTION%%%%%%%%%%%%%%%%%%%%%%%
p = configureP(pdeSolver,problem,mark,maxNrDoF,[]);
p.PDE.nu = nu2;
p.params.showIteratively = false;
p = run(p);

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = show('drawCondNr_--',saveFigures,format,p);

% % %%%%%%%%%%%%%%%%%%%%%%COMPUTE DISCRETE SOLUTION%%%%%%%%%%%%%%%%%%%%%%%
p = configureP(pdeSolver,problem,mark,maxNrDoF,[]);
p.PDE.nu = nu3;
p.params.showIteratively = false;
p = run(p);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = show('drawCondNr_-.',saveFigures,format,p);

grid off;
legend(['nu = ' sprintf('%g',nu1)],...
	['nu = ' sprintf('%g',nu2)],...
	['nu = ' sprintf('%g',nu3)]);
	


