function p = afem(p) 

mark = p.statics.mark;
refine = p.statics.refine;
solve = p.statics.solve;
createLinSys = p.statics.createLinSys;
estimate = p.statics.estimate;
enumerate = p.statics.enumerate;
postProc = p.statics.postProc;

% mark edges for refinement
p = mark(p);

% RedGreenBlue refinement of markededges
p = refine(p);

% create necessary data structures
p = enumerate(p);

% builds matrix A and RHS b
p = createLinSys(p);

% computes the solution of A*x = b
p = solve(p);

% post processes the solution x
p = postProc(p); 

% estimate error
p = estimate(p);
