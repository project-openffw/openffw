function p = drawGradU(p,lvl)

clf

pdeSolver = p.params.pdeSolver;
customDrawU = str2func([pdeSolver,'drawGradU']);

p = customDrawU(p,lvl);