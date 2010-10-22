function p = startElliptic

% problem = 'Elliptic_Square_exact';
% problem = 'Elliptic_Square_exact_JumpKappa';
% problem = 'Elliptic_HexagonalSlit_exact';
% problem = 'Elliptic_Lshape';
problem = 'Elliptic_Lshape_exact';
% problem = 'Elliptic_Waterfall_exact';

% pdeSolver = 'CR';         %'CR-Elliptic';
% pdeSolver = 'P1P0';       %'P0-P1-mixed-Elliptic';
% pdeSolver = 'P1';         %'P1-Elliptic';
% pdeSolver = 'P2';         %'P2-Elliptic';
pdeSolver = 'P3';         %'P3-Elliptic';
% pdeSolver = 'RT0P0';      %'RT0-P0-mixed-Elliptic';

maxNrDoF = 2000;

% mark = 'bulk';
% mark = 'maximum';
mark = 'uniform';
% mark = 'graded';

%%%%%%%%%%%%%%%%%%%%%COMPUTE DISCRETE SOLUTION%%%%%%%%%%%%%%%%%%%%%%%
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elliptic');
p.params.rhsIntegtrateExactDegree = 15;
p = p.statics.run(p);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(1);
% p = show('drawU',p);
% p = show('drawGrid',p);
% p = show('drawGradU',p);

figure(2);
% clf
hold all
% p.params.output.options.drawConvergenceRate = true;
p.params.errorIntegtrateExactDegree = 19;
p = show('drawError_estimatedError',p);
p = show('drawError_L2error',p);
p = show('drawError_H1semiError',p);
% hold off
% 
