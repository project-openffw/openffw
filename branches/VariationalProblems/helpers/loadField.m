function data = loadField(field,member,p,defaultValue)
% load a field from a structure
if(nargin < 4)
    defaultValue = [];
end
string1 = strcat('data=getfield(',field,',','''',member,'''',');');
string2 = 'data=defaultValue;';
eval(string1,string2);

