% Enzymatic reactions: Michaelis Menten model
% Quasi steady-state approximation

clearvars;
close all;

Vmax = 10;
Km = 0.5;

S = 0:0.01:4;
Rate1 = Vmax*S./(Km+S);

figure
hold on
plot(S,Rate1,'-b','linewidth',2);
plot(Km,Vmax/2,'or','linewidth',3);
plot([Km Km],[0 Vmax/2],'--r','linewidth',2);
plot([0 Km],[Vmax/2 Vmax/2],'--r','linewidth',2);
plot([0 S(end)],[Vmax Vmax],'--','Color',[.7 .7 .7],'linewidth',2);
plot([2*Km 2*Km],[0 2*Vmax/3],'--r','linewidth',2);
plot([0 2*Km],[2*Vmax/3 2*Vmax/3],'--r','linewidth',2);
plot(3*Km,3*Vmax/4,'or','linewidth',3);
plot([3*Km 3*Km],[0 3*Vmax/4],'--r','linewidth',2);
plot([0 3*Km],[3*Vmax/4 3*Vmax/4],'--r','linewidth',2);
plot(2*Km,2*Vmax/3,'or','linewidth',3);
axis([0 S(end) 0 12]);
set(gca,'fontsize',24);
xlabel('Substrate concentration')
ylabel('Reaction rate');
set(gca,'XTick',(0:0.5:4))
set(gca,'YTick',(0:1:12))
