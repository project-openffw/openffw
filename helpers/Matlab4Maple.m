function funcString = Matlab4Maple(func)

syms dummy real;
func = func + dummy - dummy;
func = simple(func);
funcString = char(func);

funcString = strrep(funcString,'*','.*');
funcString = strrep(funcString,'/','./');
funcString = strrep(funcString,'^','.^');

funcString = strrep(funcString,'matrix','');
funcString = strrep(funcString,'],[','];[');

funcString = strrep(funcString,',','+x-x+y-y,');
funcString = strrep(funcString,']]','+x-x+y-y]]');
funcString = strrep(funcString,'];[','+x-x+y-y];[');

funcString = strrep(funcString,'param','p.PDE.');