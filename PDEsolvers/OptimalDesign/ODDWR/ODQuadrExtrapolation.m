function QuardExtrapolation

s(1) = -0.375000589373703;
s(2) = -0.765344656478005;
s(3) = -0.796382686884807;
s(4) = -0.837568592000768;
s(5) = -0.841870191650921;
 
while length(s) >= 4
    N = length(s);
    xk_2 = s(2:N-2)-s(1:N-3);
    xk_1 = s(3:N-1)-s(1:N-3);
    xk_0 = s(4:N)  -s(1:N-3);
    s = (-xk_2.*xk_0 - xk_1.*xk_0+1).*s(2:N-2) + (1-xk_1.*xk_0).*s(3:N-1) + s(4:N);
    s
end

fprintf('Value = %.15g', s(end));