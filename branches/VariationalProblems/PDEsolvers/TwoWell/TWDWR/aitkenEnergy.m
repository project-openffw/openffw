function aitkenEnergy


s(1) = 1.17446544619082;
s(2) = 1.08013363641975;
s(3) = 1.07100610836094;
s(4) = 1.06152593292619;

while length(s) >= 3
    N = length(s);
    s = s(1:N-2)- (s(2:N-1)-s(1:N-2)).^2./(s(3:N)-2*s(2:N-1)+s(1:N-2));
end

fprintf('Value = %.15g', s(end));

% for curLvl = 1 : length(p.level)
%     t = getfield(p.level, {curLvl} , value);
%     p = setfield(p, ['level(',int2str(curLvl),').',value,'Error'], norm(t-s(end)));
% end

