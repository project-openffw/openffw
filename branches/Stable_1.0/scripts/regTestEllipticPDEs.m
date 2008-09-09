function regTestEllipticPDEs
warning off
restoredefaultpath
warning on
cd ..
addpath(fullfile(pwd,'scripts'));
addpath(fullfile(pwd,'helpers'));
addpath(fullfile(pwd,'algorithms','misc'));
maxNrDoF = 1000;



pdeSolvers = struct('z',{'RT0-P0-mixed-Elliptic','P1-Elliptic','CR-Elliptic'} );				
problems = struct('z',{'Elliptic-Exact-Template',...
					  'General-Square-exact',...
					  'Laplace-Lshape-exact',...
                      'Laplace-Mass-Square-exact'} );
                  
errTypes = struct('z',{'estimatedError','H1semiError','L2error'} );
marking = struct('z',{'bulk','uniform'});
				  
for i = 1:length(pdeSolvers)
	pdeSolver = pdeSolvers(i).z;
	for j = 1:length(problems)
        figure
		problem = problems(j).z;	
		for m = 1:length(marking)
			mark = marking(m).z;
             p = configureP(pdeSolver,problem,mark,maxNrDoF,[],'elliptic');
             p = run(p);
			for k = 1:length(errTypes)
                errType = errTypes(k).z;
               
                drawErrorString = ['drawError_',errType];
                hold all
                p = show(drawErrorString,p);
                hold off
%                 p = getConvergenceRate(p,100,errType);		
            end
        end
	end
end
		
	
