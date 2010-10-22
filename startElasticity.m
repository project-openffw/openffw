function p = startElasticity

% problem = 'Elasticity_Cooks';
% problem = 'Elasticity_Exact_Template';
problem = 'Elasticity_Lshape_exact';
% problem = 'Elasticity_Square_exact';
% problem = 'Elasticity_Square_Neumann_exact';
% problem = 'Elasticity_Template';

% pdeSolver = 'P1P1';         %'P1-P1-Elasticity';
pdeSolver = 'P1CR';       %'P1-CR-Elasticity';
% pdeSolver = 'AW';         %'Arnold-Winther-mixed-Elasticity';

maxNrDoF = 20000;

mark = 'bulk';
% mark = 'maximum';
% mark = 'uniform';
% mark = 'graded';

%%%%%%%%%%%%%%%%%%%%%%COMPUTE DISCRETE SOLUTION%%%%%%%%%%%%%%%%%%%%%%%
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elasticity');
p.PDE.nu = 0.499;
p.params.rhsIntegtrateExactDegree = 19;
p = p.statics.run(p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% p.params.output.lineWidth = 1;
% p.params.output.factor = 1000;

% figure(1);
% p = show('drawU',p);

% figure(2);
% p = show('drawU',p);
% p = show('drawGrid',p);
% p = show('drawGradU',p);

figure(2);
% p = show('drawError_estimatedError',p);
% p = show('drawError_L2error',p);
p = show('drawError_H1semiError',p);


