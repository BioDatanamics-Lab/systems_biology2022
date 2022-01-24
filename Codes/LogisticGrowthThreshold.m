% Solution of the logistic growth with a threshold equation
% Numerical method: modified Euler (Runge-Kutta, order 2)

clearvars;
close all;

% ODE
% X' = -r (1-X/T) ( 1 - X/K) X + I;

% Baseline parameters 
% r = 1;
% T = 0.4;
% K = 1;
% I = 0.0;

r = 1;
T = 0.4;
K = 1;

I = 0.02;

logthr=@(x) -r*x.*(1-x/T).*(1-x/K);


Tmax = 100;
dt = 0.01;
t = 0:dt:Tmax;

x = zeros(1,length(t));

x(1) = 1.0575;


for j=1:length(t)-1
    k1x = logthr(x(j))+I;
    ax = x(j)+k1x*dt;
    k2x = logthr(ax)+I;
    x(j+1)=x(j)+(k1x+k2x)*dt/2;
end

xx = -1:0.01:2;

figure
hold on
plot(t,x,'-b','linewidth',2);
plot([0 Tmax],[T T],'--');
plot([0 Tmax],[0 0]','--');
plot([0 Tmax],[K K],'--');
axis([0 Tmax -0.1 1.1]);
set(gca,'fontsize',20);
xlabel('t');
ylabel('V');

figure
hold on
plot([xx(1)-1 xx(end)+1],[0 0],'--','Color',[.6 .6 .6]);
plot(xx,logthr(xx)+I,'r','linewidth',2);
plot(xx,logthr(xx),'--r','linewidth',1);
axis([-0.1 1.1 -0.5 0.5]);
set(gca,'fontsize',20);
xlabel('V');
ylabel('F');

% figure(1)
% hFig = figure(1);
% set(hFig, 'Position', [40 400 1000 500]); 
% subplot(1,2,1)
% hold on
% plot(t,x,'-b','linewidth',2);
% plot([0 Tmax],[T T],'--');
% plot([0 Tmax],[0 0]','--');
% plot([0 Tmax],[K K],'--');
% axis([0 Tmax -0.1 1.1]);
% set(gca,'fontsize',20);
% xlabel('t');
% ylabel('X');
% subplot(1,2,2)
% hold on
% plot([xx(1)-1 xx(end)+1],[0 0],'--','Color',[.6 .6 .6]);
% plot(xx,logthr(xx)+I,'linewidth',2);
% axis([-0.1 1.1 -0.5 0.5]);
% set(gca,'fontsize',20);
% xlabel('X');
% ylabel('F');

% Make MVN=1 to see the phase-space diagram evolving in time ("movie").
MVN = 0;
if MVN == 1
    xx = -1:0.01:2;
    figure(101)
    k = 10;
    hFig = figure(101);
    set(hFig, 'Position', [20 100 700 500])
    for j=1:floor(length(t)/k)-1
        subplot(1,2,1)
        plot(t,x,'-b','linewidth',2);
        hold on
        plot(t(k*j),x(k*j),'or','linewidth',2);
        plot([0 Tmax],[T T],'--');
        plot([0 Tmax],[0 0]','--');
        plot([0 Tmax],[K K],'--');
        axis([0 Tmax -0.1 1.1]);
        set(gca,'fontsize',20);
        xlabel('t');
        ylabel('X');
        pause(0.0001);
        hold off  
        subplot(1,2,2)
        plot([xx(1)-1 xx(end)+1],[0 0],'--','Color',[.6 .6 .6]);
        hold on
        plot(xx,logthr(xx)+I,'linewidth',2);
        plot(x(k*j),0,'or','linewidth',2);
        axis([-0.1 1.1 -0.5 0.5]);
        set(gca,'fontsize',20);
        xlabel('X');
        pause(0.0001);
        hold off   
    end
end





