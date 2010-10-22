function drawTriangle(problem,rate1,rate2)

if strcmp(problem,'Square') %&& rate == 0.5
%     rate = 1/2;
    format rat
    str_rate = num2str(rate1);
    orient = 1;
    c = 5;
    a1 = 6 * 10^3;
    b1 = 3.5 * 10^(-3);
    x1 = [1.3*10^4,3.1*10^(-3)];
    x2 = [4.1*10^3,5.3*10^(-3)];
    maleDreieck(rate1,str_rate,a1,b1,c,...
        x1(1),x1(2),x2(1),x2(2),orient);
% elseif strcmp(problem,'Square') && rate == 1
%     rate = 1/2;
%     str_rate = num2str(rate2);
%     orient = 1;
%     c = 5;
%     a1 = 6 * 10^3;
%     b1 = 2.5 * 10^(-6);
%     x1 = [1.3*10^4,2.1*10^(-6)];
%     x2 = [5.1*10^3,4.3*10^(-6)];
%     maleDreieck(rate2,str_rate,a1,b1,c,...
%         x1(1),x1(2),x2(1),x2(2),orient);
elseif strcmp(problem,'Lshape') %&& rate == 0.5
%     rate = 1/2;
    str_rate = num2str(rate1);
    orient = 1;
    c = 5;
    a1 = 5 * 10^3;
    b1 = 2.5 * 10^(-4);
    x1 = [4*10^3,2.1*10^(-4)];
    x2 = [6*10^2,3.3*10^(-4)];
    maleDreieck(rate1,str_rate,a1,b1,c,...
        x1(1),x1(2),x2(1),x2(2),orient)
    
%     str_rate = num2str(rate2);
%     orient = 1;
%     c = 5;
%     a1 = 6 * 10^2;
%     b1 = 2.5 * 10^(-4);
%     x1 = [1.3*10^4,2.1*10^(-4)];
%     x2 = [5.1*10^3,4.3*10^(-4)];
%     maleDreieck(rate2,str_rate,a1,b1,c,...
%         x1(1),x1(2),x2(1),x2(2),orient);
end

grid off

function maleDreieck(rate,str_rate,a1,b1,c,x1,y1,x2,y2,orient)
FontSize = 10;

if orient == 1
    a2 = c * a1; b2 = (c)^(rate) * b1;
    loglog([a1;a2;a1;a1],[b1;b1;b2;b1],'k');
else
    a2 = c * a1; b2 = (c)^(-rate) * b1;
    loglog([a1;a2;a2;a1],[b1;b1;b2;b1],'k');
end

text(x1,y1,'1','FontSize',FontSize)
text(x2,y2,str_rate,'FontSize',FontSize)

h = gca;
set(h,'FontSize',FontSize)