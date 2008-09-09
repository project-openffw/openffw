addpath(fullfile(pwd,'scripts'));
addpath(fullfile(pwd,'helpers'));

% matFile = 'etaT2';

matFile = 'avg';
% matFile = 'avg_angle';
% matFile = 'uniform';
% matFile = 'uniform_SQ';
% matFile = 'etaT1';
% matFile = 'full_estimate_2_thetas';
load(matFile);


hold all;
p = show('drawError_exactH1semiError',false,'.jpg',p);
% p = show('drawError_exactL2error',false,'.jpg',p);
% p = show('drawError_estimatedError',false,'.jpg',p);
% p = show('drawGrid',false,'.jpg',p);
