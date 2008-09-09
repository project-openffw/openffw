% format = '.fig';
format = '.jpg';
% format = '.eps';

s = '_';

% problem = 'Elasticity-Lshape';
% problem = 'Elasticity-Square';
problem = 'Elasticity-Cooks';

% maxNrDoF = 2000;
maxNrDoF = 500;
DoFString = sprintf('%d',maxNrDoF);

% mark = 'bulk'	
mark = 'uniform'

% errType = 'drawError_exactL2error';
% errType = 'drawError_exactH1semiError';
errType = 'drawError_estimatedError';

% nu1 = 0.3;
% nu2 = 0.49;
% nu3 = 0.499;


nu1 = 0.3;
nu2 = 0.4999;
nu3 = 0.49999999;



% nu1 = 0.4999;
% nu2 = 0.499999;
% nu3 = 0.4999999999;

exactH1semiErrorToleranceNU1 = 1;
exactH1semiErrorToleranceNU2 = 10;
exactH1semiErrorToleranceNU3 = 100;
% exactH1semiErrorToleranceNU1 = 1e-7;
% exactH1semiErrorToleranceNU2 = 1e-7;
% exactH1semiErrorToleranceNU3 = 1e-7;

saveFigures = false;

clf;

pdeSolver = 'Arnold-Winther-mixed-Elasticity';
LockingComparison;
saveas(gcf,['LockingComparison' s pdeSolver s problem s mark s DoFString format]);
clf

% 
% pdeSolver = 'P1-CR-Elasticity';
% LockingComparison;
% saveas(gcf,['LockingComparison' s pdeSolver s problem s mark s DoFString format]);
% clf
% 
% pdeSolver = 'P1-P1-Elasticity';
% % exactH1semiErrorToleranceNU1 = 10;
% % exactH1semiErrorToleranceNU2 = 1000;
% % exactH1semiErrorToleranceNU3 = 10000000;
% 
% exactH1semiErrorToleranceNU1 = 1e-7;
% exactH1semiErrorToleranceNU2 = 1e-5;
% exactH1semiErrorToleranceNU3 = 1e-3;
% LockingComparison;
% saveas(gcf,['LockingComparison' s pdeSolver s problem s mark s DoFString format]);
% clf
% 
% 



% 
% mark = 'bulk'	
% % mark = 'uniform'
% 
% nu1 = 0.3;
% nu2 = 0.49;
% nu3 = 0.499;
% 
% exactH1semiErrorToleranceNU1 = 1;
% exactH1semiErrorToleranceNU2 = 1;
% exactH1semiErrorToleranceNU3 = 1;
% 
% saveFigures = false;
% 
% clf;
% 
% pdeSolver = 'Arnold-Winther-mixed-Elasticity-Paper';
% LockingComparison;
% saveas(gcf,['LockingComparison' s pdeSolver s problem s mark s DoFString format]);
% clf
% 
% 
% pdeSolver = 'P1-CR-Elasticity';
% LockingComparison;
% saveas(gcf,['LockingComparison' s pdeSolver s problem s mark s DoFString format]);
% clf
% 
% pdeSolver = 'P1-P1-Elasticity';
% exactH1semiErrorToleranceNU2 = 10;
% exactH1semiErrorToleranceNU3 = 100;
% LockingComparison;
% saveas(gcf,['LockingComparison' s pdeSolver s problem s mark s DoFString format]);
% clf





