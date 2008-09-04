% format = '.fig';
% format = '.jpg';
format = '.eps';

s = '_';

% problem = 'Elasticity-Lshape';
problem = 'Elasticity-Square';

% maxNrDoF = 20000;
maxNrDoF = 501;
DoFString = sprintf('%d',maxNrDoF);

% mark = 'bulk'	
mark = 'uniform'

nu = 0.3;

exactH1semiErrorTolerance = 1;
exactL2errorTolerance = 1e-8;

saveFigures = false;

% errorType = 'estimatedError';
% clf;
% globalMark = mark;
% mark = 'uniform'
% ErrorComparison;
% mark = 'bulk'	
% ErrorComparison;
% saveas(gcf,[errorType s problem s mark s DoFString format]);
% clf;
% mark = globalMark;

errorType = 'exactL2error';
ErrorComparison;
saveas(gcf,[errorType s problem s mark s DoFString format]);
clf

% errorType = 'exactH1semiError';
% ErrorComparison;
% saveas(gcf,[errorType s problem s mark s DoFString format]);
% clf



mark = 'bulk'	
% mark = 'uniform'

nu = 0.3;

exactH1semiErrorTolerance = 1;
exactL2errorTolerance = 1e-6;

saveFigures = false;


errorType = 'exactL2error';
ErrorComparison;
saveas(gcf,[errorType s problem s mark s DoFString format]);
clf

% errorType = 'exactH1semiError';
% ErrorComparison;
% saveas(gcf,[errorType s problem s mark s DoFString format]);
% clf




