
warning off
restoredefaultpath
warning on
addpath(fullfile(pwd,'scripts'));
addpath(fullfile(pwd,'helpers'));

hold all;


pdeSolver = 'Arnold-Winther-mixed-Elasticity';
% % %%%%%%%%%%%%%%%%%%%%%%COMPUTE DISCRETE SOLUTION%%%%%%%%%%%%%%%%%%%%%%%
p = configureP(pdeSolver,problem,mark,maxNrDoF,[]);
p.PDE.nu = nu;
p.params.exactL2errorTolerance = exactL2errorTolerance;
p.params.exactH1semiErrorTolerance = exactH1semiErrorTolerance;
p.params.showIteratively = false;
p = run(p);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = show(['drawError_' errorType '_NoSlope_-'],saveFigures,format,p);

pdeSolver = 'P1-P1-Elasticity';
% % %%%%%%%%%%%%%%%%%%%%%%COMPUTE DISCRETE SOLUTION%%%%%%%%%%%%%%%%%%%%%%%
p = configureP(pdeSolver,problem,mark,maxNrDoF,[]);
p.PDE.nu = nu;
p.params.exactL2errorTolerance = exactL2errorTolerance;
p.params.exactH1semiErrorTolerance = exactH1semiErrorTolerance;
p.params.showIteratively = false;
p = run(p);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = show(['drawError_' errorType '_NoSlope_--'],saveFigures,format,p);

pdeSolver = 'P1-CR-Elasticity';
% % %%%%%%%%%%%%%%%%%%%%%%COMPUTE DISCRETE SOLUTION%%%%%%%%%%%%%%%%%%%%%%%
p = configureP(pdeSolver,problem,mark,maxNrDoF,[]);
p.PDE.nu = nu;
p.params.exactL2errorTolerance = exactL2errorTolerance;
p.params.exactH1semiErrorTolerance = exactH1semiErrorTolerance;
p.params.showIteratively = false;
p = run(p);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = show(['drawError_' errorType '_NoSlope_-.'],saveFigures,format,p);

grid off;
legend('AW','P1','KS')


