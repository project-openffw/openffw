function p = computeSolution(p)

solver = str2func(p.params.solver);
p = solver(p);
