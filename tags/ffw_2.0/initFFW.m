function p = initFFW(pdeSolver,problem,mark,maxNrDoF,problemType,solver,refine,estimate,trafo)
% Initialisation of the FFW.
% Creates the structure p with all parameters and functions handles.
% All parameters are saved in
%         p.params
% while the function handles are stored in
%         p.statics
%
% input:
% pdeSolver   - Type of finite element method ('P1','P2','P3','CR',...)
% problem     - Definition of a problem (string of the name of a file in
%               ./problem/elliptic or ./problem/elasticity)
% mark        - marking strategy ('uniform','bulk','maximum','graded')
% maxNrDoF    - integer indicating when to stop ( maximum number of degrees
%               of freedom)
% problemType - Type of the problem ('elliptic' or 'elasticity')
% solver      - optional [{'direct'},'multigrid']
% refine      - optional [{'redGreenBlue'}] for further developement
% estimate    - optional [{'estimate'}] for further developement
% trafo       - optional [{'biaffin','isoparametric'}] for further developement

% Copyright 2007 Jan Reininghaus, David Guenther, 
%                Andreas Byfut, Joscha Gedicke
%
% This file is part of FFW.
%
% FFW is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% FFW is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

if nargin < 6
    solver = 'direct';
    refine = 'redGreenBlue';
    estimate = 'estimate';
    trafo = 'affin';
elseif nargin < 7
    refine = 'redGreenBlue';
    estimate = 'estimate';
    trafo = 'affin';
elseif nargin < 8
    estimate = 'estimate';
    trafo = 'affin';
elseif nargin < 9
    trafo = 'affin';
end
    

%% add paths 
warning off;
restoredefaultpath;
warning on;
curPath = pwd;
subDir = [findstr('algorithms',curPath),findstr('evaluation',curPath),...
          findstr('helpers',curPath),findstr('PDEsolvers',curPath),...
          findstr('problems',curPath),findstr('results',curPath),findstr('scripts',curPath)];
if ~isempty(subDir)
    curPath = curPath(1:subDir-2);
end
addpath(genpath(curPath))

%% show PDEsolver used
disp(pdeSolver)

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

%% load general algorithms as specified
p.statics.mark = str2func(p.params.mark);
p.statics.refine = str2func(p.params.refine);
p.statics.solve = str2func(p.params.solver);

p.statics.closure = str2func('closure');

%% load generic enumerate functions
p.statics.enum = [];
p.statics.enum.getE4n = str2func('getE4n');
p.statics.enum.getEd4n = str2func('getEd4n');
p.statics.enum.getN4ed = str2func('getN4ed');
p.statics.enum.getEd4e = str2func('getEd4e');
p.statics.enum.getE4ed = str2func('getE4ed');
p.statics.enum.getDbEdges = str2func('getDbEdges');
p.statics.enum.getNbEdges = str2func('getNbEdges');
p.statics.enum.getArea4e = str2func('getArea4e');
p.statics.enum.getArea4n = str2func('getArea4n');
p.statics.enum.getLength4ed = str2func('getLength4ed');
p.statics.enum.getMidPoint4e = str2func('getMidPoint4e');
p.statics.enum.getMidPoint4ed = str2func('getMidPoint4ed');
p.statics.enum.getNormals4e = str2func('getNormals4e');
p.statics.enum.getTangents4e = str2func('getTangents4e');
p.statics.enum.getAngles4e = str2func('getAngles4e');
p.statics.enum.getAngle4n = str2func('getAngle4n');
p.statics.enum.getNormals4NbEd = str2func('getNormals4NbEd');
p.statics.enum.getNormals4DbEd = str2func('getNormals4DbEd');
p.statics.enum.getNormals4ed = str2func('getNormals4ed');

%% specify transformation formula to be used
eval(['p = transformation_',trafo,'(p);']);

%% load PDEsolver specific functions
estimate = p.params.estimate;
if(~exist([pdeSolver,estimate,'.m'],'file'))
    p.statics.estimate = @STUBestimate;
    if(~strcmp(mark,'uniform') && ~strcmp(mark,'graded'))
        warning('Adaptive refinement needs a valid error estimater');
    end
else
    p.statics.estimate = str2func([pdeSolver,estimate]);
end
p.statics.createLinSys = str2func([pdeSolver,'createLinSys']);
p.statics.enumerate = str2func([pdeSolver,'enumerate']);
p.statics.postProc = str2func([pdeSolver,'postProc']);
p.statics.init = str2func([pdeSolver,'init']);
  
%% load problem data
getProblem = str2func(p.params.problem.name);
p = getProblem(p);

%% load geometry data
curGeom = p.problem.geom;

% Check if geometry exists
pathPrefix = fullfile(curPath,'problems','geometries',curGeom);
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

% Element-Type Specific Initialization 
p = p.statics.init(p);

%% default integration parameters
p.params.integrationDegrees.exactError = 9;
% Further FEM-based integration parameters are set in the FEM specific
% initialization.







