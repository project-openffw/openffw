function p = computeSolution(p)
% Main computation loop. Computes the AFEM loop until stopping criterion is
% reached and specifies functions to evaluate the error

% Set initial values
p.level(1).level = 1;
p.level(1).estimatedError = 0;

% load abort criteria
maxLevel = loadField('p.params','maxLevel',p,1000);
maxNrDoF = p.params.maxNrDoF;

curLevel = 1;
p.level(1).level = curLevel;

% refine first levels
refineFirstLevel = p.params.modules.mark.refineFirstLevel;
mark = p.statics.mark;
refine = p.statics.refine;
enumerate = p.statics.enumerate;
while( curLevel < refineFirstLevel )
    p = mark(p);
    p = refine(p);
    p = enumerate(p);
    curLevel = curLevel+1;
    p.level(curLevel).level = curLevel;
end

nrDoF = p.level(end).nrDoF;

fprintf('\n');
while(nrDoF < maxNrDoF && curLevel <= maxLevel)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = afem(p);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nrDoF = p.level(end).nrDoF;
    estimatedError = p.level(end).estimatedError;
    fprintf('Level = %d \t Error = %.2g \t DoF = %d\n',curLevel+1,estimatedError,nrDoF);
    curLevel = curLevel + 1;
    p.level(end).level = curLevel;
    maxNrDoF = p.params.maxNrDoF;
end

