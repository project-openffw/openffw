
mu1 = 1;
mu2 = 2;

lambda = 0.1;
xi = 0.5;

theta = 0;

t1 = sqrt(2*lambda*(mu1/mu2));
t2 = (mu2/mu1)*t1;

mu2*t1-mu1*t2

g1 = @(t) 0.5*mu2*t^2+lambda*xi*(mu1-mu2);
g2 = @(t) t1*mu2*(t-0.5*t1)+lambda*xi*(mu1-mu2);
g3 = @(t) 0.5*mu1*t^2 + 0.5*mu1*t2^2 - 0.5*t1^2*mu2+lambda*xi*(mu1-mu2);

h1 = @(t) t^2/(2*mu2) - lambda*xi*(mu1-mu2);
h2 = @(t) t^2/(2*mu1) -mu1*t2^2/2 + mu2*t1^2/2 - lambda*xi*(mu1-mu2);

f1 = @(t) 0.5*t^2/(mu2+theta) - lambda*xi*(mu1-mu2);
f2 = @(t) 1/(2*theta)*(mu2*t1-t)^2 + mu2*t1^2/2  - lambda*xi*(mu1-mu2);
f3 = @(t) 0.5*t^2/(mu1+theta) - mu1*t2^2/2 + mu2*t1^2/2 - lambda*xi*(mu1-mu2);

dg1 = @(t) mu2*t;
dg2 = @(t) t1*mu2*ones(length(t),1);
dg3 = @(t) mu1*t;

dh1 = @(t) t/(mu2);
dh2 = @(t) t/(mu1);

xMin = 0;
xMax = 1.5;
yMin = 0;
yMax = 1;

fontSize = 18;
markerSize = 10;

energyDensity = false;
energyDensityConj = false;
energyDensityConjReg = false;
derEnergyDensity = true;
derEnergyDensityConj = false;

if energyDensity
    figure(1)
    clf
    hold all
    ezplot(g1,[xMin xMax,yMin yMax])
    ezplot(g2,[xMin xMax,yMin yMax])
    ezplot(g3,[xMin xMax,yMin yMax])
    line([t1,t1],[yMin,yMax],'Linestyle','--','Color','k')
    plot(t1,yMin,'Marker','x','Markersize',markerSize,'Color','k')
    text(t1,yMin-0.05,'t_1','Fontsize',fontSize)
    line([t2,t2],[yMin,yMax],'Linestyle','--','Color','k')
    plot(t2,yMin,'Marker','x','Markersize',markerSize,'Color','k')
    text(t2,yMin-0.05,'t_2','Fontsize',fontSize)
    legend('0.5 \mu_2 t^2+\lambda \xi (\mu_1-\mu_2)','t_1 \mu_2 (t-0.5 t_1)+\lambda \xi (\mu_1-\mu_2)','0.5 \mu_1 t^2 + 0.5 \mu_1 t_2^2 - 0.5 t_1^2 \mu_2+\lambda \xi (\mu_1-\mu_2)','Location','Northwest')
    title('')
    xlabel('')
    set(gca,'Fontsize',fontSize)
end

if energyDensityConj
    figure(2)
    clf
    hold all
    h = ezplot(h1,[xMin xMax,yMin yMax])
    set(h,'LineWidth',3,'Linestyle','-');
    h = ezplot(h2,[xMin xMax,yMin yMax])
    set(h,'LineWidth',3,'Linestyle','--');
    line([t1*mu2,t1*mu2],[yMin,yMax],'Linestyle','--','Color','k')
    plot(t1*mu2,yMin,'Marker','x','Markersize',markerSize,'Color','k')
    text(t1*mu2,yMin-0.05,'t_1\mu_2','Fontsize',fontSize)
    legend('t^2/(2 \mu_2) + \lambda \xi (\mu_1-\mu_2)','t^2/(2 \mu_1) -\mu_1 t_2^2/2 + \mu_2 t_1^2/2 + \lambda \xi (\mu_1-\mu_2)','Location','Southeast')
    title('')
    xlabel('')
    set(gca,'Fontsize',fontSize)
end

if energyDensityConjReg
    figure(3)
    clf
    hold all
    h = ezplot(f1,[xMin xMax,yMin yMax]);
    set(h,'LineWidth',3,'Linestyle','-');
    h = ezplot(f2,[xMin xMax,yMin yMax]);
    set(h,'LineWidth',3,'Linestyle','--');
    h = ezplot(f3,[xMin xMax,yMin yMax]);
    set(h,'LineWidth',3,'Linestyle',':');
    line([t1*(mu2+theta),t1*(mu2+theta)],[yMin,yMax],'Linestyle','--','Color','k')
    plot(t1*(mu2+theta),yMin,'Marker','x','Markersize',markerSize,'Color','k')
    text(t1*(mu2+theta)-0.15,yMin-0.1,'t_1(\mu_2+\theta)','Fontsize',fontSize)
    line([t2*(mu1+theta),t2*(mu1+theta)],[yMin,yMax],'Linestyle','--','Color','k')
    plot(t2*(mu1+theta),yMin,'Marker','x','Markersize',markerSize,'Color','k')
    text(t2*(mu1+theta)+0.05,yMin-0.1,'t_2(\mu_1+\theta)','Fontsize',fontSize)
    legend('0.5 t^2/(\mu_2+\theta) - \lambda \xi (\mu_1-\mu_2)','1/(2 \theta) (\mu_2 t_1-t)^2 + \mu_2 t_1^2/2  - \lambda \xi (\mu_1-\mu_2)','0.5 t^2/(\mu_1+\theta) - \mu_1 t_2^2/2 + \mu_2 t_1^2/2 - \lambda \xi (\mu_1-\mu_2)','Location','Southeast')
    title('')
    xlabel('')
    set(gca,'Fontsize',fontSize)
end

if derEnergyDensity
    figure(4)
    clf
    hold all
    points1 = xMin:1/100:t1;
    evalDg1 = dg1(points1);
%     ezplot(dg1,[xMin t1,yMin yMax])
    h = plot(points1,evalDg1);
    set(h,'LineWidth',3,'Linestyle','-');
    points2 = t1:1/100:t2;
    evalDg2 = dg2(points2);
%     ezplot(dg2,[t1 t2,yMin yMax])
    h = plot(points2,evalDg2);
    set(h,'LineWidth',3,'Linestyle','-');
    points3 = t2:1/100:xMax;
    evalDg3 = dg3(points3);
%     ezplot(dg3,[t2 xMax,yMin yMax])
    h = plot(points3,evalDg3);
    set(h,'LineWidth',3,'Linestyle','-');
    line([t1,t1],[yMin,yMax],'Linestyle','--','Color','k')
    plot(t1,yMin,'Marker','x','Markersize',markerSize,'Color','k')
    text(t1,yMin-0.05,'t_1','Fontsize',fontSize)
    line([t2,t2],[yMin,yMax],'Linestyle','--','Color','k')
    plot(t2,yMin,'Marker','x','Markersize',markerSize,'Color','k')
    text(t2,yMin-0.05,'t_2','Fontsize',fontSize)
    legend('\mu_2 t','t_1 \mu_2','\mu_1 t','Location','Northwest')
    title('')
    xlabel('')
    set(gca,'Fontsize',fontSize) 
    axis auto
end

if derEnergyDensityConj
       figure(5)
    clf
    hold all
    points1 = xMin:1/100:mu2*t1;
    evalDh1 = dh1(points1);
%     ezplot(dg1,[xMin t1,yMin yMax])
    h = plot(points1,evalDh1);
    set(h,'LineWidth',3,'Linestyle','-');
    points2 = mu2*t1:1/100:xMax;
    evalDh2 = dh2(points2);
%     ezplot(dg2,[t1 t2,yMin yMax])
    h = plot(points2,evalDh2);
    set(h,'LineWidth',3,'Linestyle','-');
    line([t1*mu2,t1*mu2],[yMin,yMax],'Linestyle','--','Color','k')
    plot(t1*mu2,yMin,'Marker','x','Markersize',markerSize,'Color','k')
    text(t1*mu2,yMin-0.05,'t_1\mu_2','Fontsize',fontSize)
    legend('t/\mu_2','t/\mu_1','Location','Southeast')
    title('')
    xlabel('')
    set(gca,'Fontsize',fontSize) 
    axis auto
end