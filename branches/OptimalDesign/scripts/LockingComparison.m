
warning off
restoredefaultpath
warning on
addpath(fullfile(pwd,'scripts'));
addpath(fullfile(pwd,'helpers'));

hold all;

% % %%%%%%%%%%%%%%%%%%%%%%COMPUTE DISCRETE SOLUTION%%%%%%%%%%%%%%%%%%%%%%%
p = configureP(pdeSolver,problem,mark,maxNrDoF,[]);
p.PDE.nu = nu1;
p.params.exactH1semiErrorTolerance = exactH1semiErrorToleranceNU1;
p.params.exactL2errorTolerance = exactH1semiErrorToleranceNU1;
p.params.showIteratively = false;
p = run(p);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
errTypeEXt = [errType '_NoSlope_-'];
p = show(errTypeEXt,saveFigures,format,p);


% % %%%%%%%%%%%%%%%%%%%%%%COMPUTE DISCRETE SOLUTION%%%%%%%%%%%%%%%%%%%%%%%
p = configureP(pdeSolver,problem,mark,maxNrDoF,[]);
p.PDE.nu = nu2;
p.params.exactH1semiErrorTolerance = exactH1semiErrorToleranceNU2;
p.params.exactL2errorTolerance = exactH1semiErrorToleranceNU2;
p.params.showIteratively = false;
p = run(p);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
errTypeEXt = [errType '_NoSlope_--'];
p = show(errTypeEXt,saveFigures,format,p);

% % %%%%%%%%%%%%%%%%%%%%%%COMPUTE DISCRETE SOLUTION%%%%%%%%%%%%%%%%%%%%%%%
p = configureP(pdeSolver,problem,mark,maxNrDoF,[]);
p.PDE.nu = nu3;
p.params.exactH1semiErrorTolerance = exactH1semiErrorToleranceNU3;
p.params.exactL2errorTolerance = exactH1semiErrorToleranceNU3;
p.params.showIteratively = false;
p = run(p);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
errTypeEXt = [errType '_NoSlope_-.'];
p = show(errTypeEXt,saveFigures,format,p);

grid off;
legend(['nu = ' sprintf('%g',nu1)],...
	['nu = ' sprintf('%g',nu2)],...
	['nu = ' sprintf('%g',nu3)]);


