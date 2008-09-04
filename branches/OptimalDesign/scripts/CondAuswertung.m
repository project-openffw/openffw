% format = '.fig';
% format = '.jpg';
format = '.eps';

s = '_';

problem = 'Elasticity-Lshape';
% problem = 'Elasticity-Square';

% maxNrDoF = 12000;
maxNrDoF = 1400;
DoFString = sprintf('%d',maxNrDoF);

% mark = 'bulk'	
mark = 'uniform'

nu1 = 0.3;
nu2 = 0.49;
nu3 = 0.4999;

saveFigures = false;

clf;

% pdeSolver = 'Arnold-Winther-mixed-Elasticity-Paper';
% CondComparison;
% saveas(gcf,['CondComparison' s pdeSolver s problem s mark s DoFString format]);
% clf
% 
% % 
% pdeSolver = 'P1-CR-Elasticity';
% CondComparison;
% saveas(gcf,['CondComparison' s pdeSolver s problem s mark s DoFString format]);
% clf

pdeSolver = 'P1-P1-Elasticity';
CondComparison;
saveas(gcf,['CondComparison' s pdeSolver s problem s mark s DoFString format]);
clf



mark = 'bulk'	
% mark = 'uniform'


clf;

pdeSolver = 'Arnold-Winther-mixed-Elasticity-Paper';
CondComparison;
saveas(gcf,['CondComparison' s pdeSolver s problem s mark s DoFString format]);
clf


pdeSolver = 'P1-CR-Elasticity';
CondComparison;
saveas(gcf,['CondComparison' s pdeSolver s problem s mark s DoFString format]);
clf

pdeSolver = 'P1-P1-Elasticity';
CondComparison;
saveas(gcf,['CondComparison' s pdeSolver s problem s mark s DoFString format]);
clf





