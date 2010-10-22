function p = initFFW(pdeSolver,problem,mark,maxNrDoF,problemType,solver,refine,estimate,errorOutput)

if(nargin < 6)
    solver = 'direct';
end
if(nargin < 7)
    refine = 'redGreenBlue';
end
if(nargin < 8)
    estimate = 'estimate';
end
if(nargin < 9)
    errorOutput = true;
end

%% add paths 
warning off;
restoredefaultpath;
warning on;
addpath(genpath(pwd))

%% show PDEsolver used
if errorOutput
    disp(pdeSolver)
end

%% initialize p with defaultParameters and values given
p = [];   
p.params.mark = mark;
p.params.problem.name = problem;
p.params.pdeSolver = pdeSolver;
p.params.maxNrDoF = maxNrDoF;
p.params.problem.type = problemType;
p.params.solver = solver;
p.params.refine = refine;
p.params.estimate = estimate;
p.params.errorOutput = errorOutput;

% mark generics
p.params.modules.mark.refineFirstLevel = 1;

%% load general algorithms as specified
p.statics.mark = str2func(p.params.mark);
p.statics.refine = str2func(p.params.refine);
p.statics.solve = str2func(p.params.solver);

%% load PDEsolver specific functions
estimate = p.params.estimate;
if(~exist([pdeSolver,estimate,'.m'],'file'))
    p.statics.estimate = @STUBestimate;
    if(~strcmp(mark,'uniform') && ~strcmp(mark,'graded'))
        error('Adaptive refinement needs a valid error estimater');
    end
else
    p.statics.estimate = str2func([pdeSolver,estimate]);
end
p.statics.createLinSys = str2func([pdeSolver,'createLinSys']);
p.statics.enumerate = str2func([pdeSolver,'enumerate']);
p.statics.postProc = str2func([pdeSolver,'postProc']);
p.statics.run = str2func([pdeSolver,'run']);
  
%% load problem data
getProblem = str2func(p.params.problem.name);
p = getProblem(p);

%% load geometry data
curGeom = p.problem.geom;

% Check if geometry exists
pathPrefix = fullfile(pwd,'problems','geometries',curGeom);
if(exist(pathPrefix,'dir' ) ~= 7)
		error('geometry does not exist! (spelling?)');
end

% Load geometry
prefix = [curGeom,'_'];
p.level(1).geom.c4n = load(fullfile(pathPrefix,[prefix,'c4n.dat']));
p.level(1).geom.n4e = load(fullfile(pathPrefix,[prefix,'n4e.dat']));
p.level(1).geom.Db  = load(fullfile(pathPrefix,[prefix,'Db.dat']));
eval('p.level(1).geom.Nb=load(fullfile(pathPrefix,[prefix,''Nb.dat'']));',...
     'p.level(1).geom.Nb=[];');
 
% Load specific problemType_init
initPDEtypeFunc = str2func([problemType,'_init']);
p = initPDEtypeFunc(p);








