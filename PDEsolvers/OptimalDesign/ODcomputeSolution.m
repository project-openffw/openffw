function p = ODcomputeSolution(p)

solver = str2func(p.params.solver);
p = solver(p);
