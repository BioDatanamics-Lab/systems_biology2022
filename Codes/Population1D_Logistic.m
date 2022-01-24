% 2022-01-24

% Population increase: 1D quadratic dynamics - Logistic equation

clearvars;
close all;

% Parameters


a = 0.21;
c = 0.2;
N = 100;

% Time definitions

Tmax = 1000;
dt = 0.01;
t = 0:dt:Tmax;

% Numerical Simulations

I = zeros(1,length(t));

I(1) = 2;

for j=1:length(t)-1
    k1i = (a-c)*I(j)-a*I(j)^2/N;
    ai = I(j)+k1i*dt;
    k2i = (a-c)*ai-a*ai^2/N;
    I(j+1) = I(j)+(k1i+k2i)*dt/2;
end

figure
hold on
plot(t,I,'-b','linewidth',2);
axis([0 Tmax 0 10]);
set(gca,'fontsize',24);
xlabel('t   [au]');
ylabel('I');

xi = 0:0.01:10;
F = (a-c)*xi-a*xi.^2/N;

figure
hold on
plot(xi,F,'-b','linewidth',2);
plot(xi,zeros(1,length(xi)),'--k','linewidth',1);
axis([0 10 -0.01 0.02]);
set(gca,'fontsize',24);
xlabel('I');
ylabel('F');

