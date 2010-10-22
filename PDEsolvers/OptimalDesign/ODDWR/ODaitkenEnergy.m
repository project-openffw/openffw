function aitkenEnergy

s(1) = 1.04402146073422;
s(2) = 1.07726625068036;
s(3) = 1.09142992161587;
s(4) = 1.0934442873401;
%s=d;

while length(s) >= 3
    N = length(s);
    s = s(1:N-2)- (s(2:N-1)-s(1:N-2)).^2./(s(3:N)-2*s(2:N-1)+s(1:N-2));
    %s
end

fprintf('Value = %.15g', s(end));

% for curLvl = 1 : length(p.level)
%     t = getfield(p.level, {curLvl} , value);
%     p = setfield(p, ['level(',int2str(curLvl),').',value,'Error'], norm(t-s(end)));
% end

