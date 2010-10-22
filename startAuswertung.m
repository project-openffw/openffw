function p = startAuswertung

%% visualizing/post processing of the solution

 hold all

% load results/Test7/TwoSquare_ODP2_h-dep_estimate_4000_bulk;
% disp(['TwoSquare_ODP2_h-dep_estimate_4000_bulk'])
 load results/Test7/TwoSquare_ODP2_h-dep_estimate_4000_uniform;
 disp(['TwoSquare_ODP2_h-dep_estimate_4000_uniform'])


 p.params.output.saveFigures = true

 figure(1);
 p.params.output.name = 'estimated';
 p = show('drawError_estimatedError',p); 
 figure(2);
 p.params.output.name = 'L2';
 p = show('drawError_L2errorDisplacement',p);
 figure(3);
 p.params.output.name = 'H1semi';
 p = show('drawError_L2errorGradU',p);
 figure(4);
 p.params.output.name = 'NonLinear';
 p = show('drawError_L2errorPhminusP0',p);
 figure(5);
 p.params.output.name = 'Energy';
 p = show('drawError_errorEnergy',p);
