% 2022-01-24

% Population increase: 1D linear dynamics
% 
% The population S can be interpreted as the infected population I

clearvars;
close all;

% Parameters

Delta = 0.1;
a = 0.6;
mu = 0.5;

% Time definitions

Tmax = 100;
dt = 0.01;
t = 0:dt:Tmax;

% Numerical Simulations

S = zeros(1,length(t));

S(1) = 0.1;

for j=1:length(t)-1
    k1s = Delta+(a-mu)*S(j);
    as = S(j)+k1s*dt;
    k2s = Delta+(a-mu)*as;
    S(j+1) = S(j)+(k1s+k2s)*dt/2;
end

figure
hold on
plot(t,S,'-b','linewidth',2);
axis([0 Tmax 0 1.2]);
set(gca,'fontsize',24);
xlabel('t   [au]');
ylabel('S');