function p = conditionNumber(p)


level = p.level(end).level;
refineFirstLevel = p.params.modules.mark.refineFirstLevel;

for curLvl = refineFirstLevel+1:level
    p.level(curLvl).conditionNr = condest(p.level(curLvl).A);
end
    
    
