% 2022-02-03
% Chemical Reaction Networks
% Ingalls BP, Mathematical models in systems biology, MIT Press, 2013

% Numerical method: Modified Euler (Runge-Kutta order 2, used even for
% linear equations for uniformity purposes)

clearvars;
close all;

lightblueish = [.4 .6 .9];
lightcoral = [0.94 0.5 0.5];
lightsalmon = [0.9 0.5 0.4];
lightgray = [.7 .7 .7];
darkgray = [.3 .3 .3];
darkgray2 = [.1 .1 .1];
mediumacquamarine = [0.4 0.8 0.6];

% [t] = sec

MDL = 5;
        % 1: Decay (Example I, eq. 2.2)
        % 2: Production & Decay (Example II, eq. 2.5)
        % 3: Irreversible converstion (Example III, eq. )
        % 4: Reversible conversion (Example IV, eqs. 2.8-2.9)
        % 5: Network example (eq. 2.17)
        %    ChemicalNetwork02_Ingalls.m
        % 6: Network example (eq. 2.19) - Rapid equilibrium assumption
        %    ChemicalNetwork02_Ingalls.m

if MDL == 1

    % Parameters
    
    k = 1;
    
    % Time definitions
    
    Tmax = 20;
    dt = 0.01;
    t = 0:dt:Tmax;
    
    % Numerical solution
    
    a = zeros(1,length(t));
    
    a(1) = 1;
    
    for j=1:length(t)-1
        kaux1a = -k*a(j);
        auxa = a(j)+kaux1a*dt;
        kaux2a = -k*auxa;
        a(j+1) = a(j)+(kaux1a+kaux2a)*dt/2;
    end
   
    
    figure
    hold on
    plot(t,a,'-b','linewidth',2);
    plot([0 Tmax],[0 0],'-k');
    plot(t,a,'-b','linewidth',2);
    axis([0 Tmax -0.2 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [sec]');
    ylabel('Conc  [mM]');
    legend('a');
    
elseif MDL == 2
    
    % Parameters
    
    k0 = 5;
    k1 = 10;
    
    % Time definitions
    
    Tmax = 20;
    dt = 0.01;
    t = 0:dt:Tmax;
    
    % Numerical solution
    
    a = zeros(1,length(t));
    
    a(1) = 1;
    
    for j=1:length(t)-1
        kaux1a = k0-k1*a(j);
        auxa = a(j)+kaux1a*dt;
        kaux2a = k0-k1*auxa;
        a(j+1) = a(j)+(kaux1a+kaux2a)*dt/2;
    end
 
    figure
    hold on
    plot(t,a,'-b','linewidth',2);
    plot([0 Tmax],[0 0],'-k');
    plot(t,a,'-b','linewidth',2);
    axis([0 Tmax -0.2 2.5]);
    set(gca,'fontsize',24);
    xlabel('t  [sec]');
    ylabel('Conc  [mM]');
    legend('a');
    
elseif MDL == 3
    
    % Parameters
    
    k = 10;
    
    % Time definitions
    
    Tmax = 20;
    dt = 0.01;
    t = 0:dt:Tmax;
    
    % Numerical solution
    
    a = zeros(1,length(t));
    
    a(1) = 1;
    b(1) = 0;
    
    for j=1:length(t)-1
        kaux1a = -k*a(j);
        auxa = a(j)+kaux1a*dt;
        kaux2a = -k*auxa;
        a(j+1) = a(j)+(kaux1a+kaux2a)*dt/2;
    end
   
    b = a(1)+b(1)-a;
    
    figure
    hold on
    plot(t,a,'-b','linewidth',2);
    plot(t,b,'-r','linewidth',2);
    plot([0 Tmax],[0 0],'-k');
    plot(t,a,'-b','linewidth',2);
    axis([0 Tmax -0.2 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [sec]');
    ylabel('Conc  [mM]');
    legend('a','b');
    
elseif MDL == 4
    
        % Parameters
    
    kp = 1;
    km = 10;
    
    % Time definitions
    
    Tmax = 10;
    dt = 0.01;
    t = 0:dt:Tmax;
    
    % Numerical solution
    
    a = zeros(1,length(t));
    b = zeros(1,length(t));
    
    a(1) = 1;
    b(1) = 0;
    
    for j=1:length(t)-1
        kaux1a = km*b(j)-kp*a(j);
        kaux1b = kp*a(j)-km*b(j);
        auxa = a(j)+kaux1a*dt;
        auxb = b(j)+kaux1b*dt;
        kaux2a = km*auxb-kp*auxa;
        kaux2b = kp*auxa-km*auxb;        
        a(j+1) = a(j)+(kaux1a+kaux2a)*dt/2;
        b(j+1) = b(j)+(kaux1b+kaux2b)*dt/2;
    end
   
    
    figure
    hold on
    plot(t,a,'-b','linewidth',2);
    plot(t,b,'-r','linewidth',2);
    plot([0 Tmax],[0 0],'-k');
    plot(t,a,'-b','linewidth',2);
    axis([0 Tmax -0.2 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [sec]');
    ylabel('Conc  [mM]');
    legend('a','b');
    
elseif MDL == 5
    
    % Parameters
% 
%     k1 = 3;
%     k2 = 2;
%     k3 = 2.5;
%     k4 = 3;
%     k5 = 4;
    
    k1 = 1;
    k2 = 0.1;
    k3 = 1;
    k4 = 2;
    k5 = 3;

    % Time definitions

    Tmax = 100;
    dt = 0.01;
    t = 0:dt:Tmax;

    a = zeros(1,length(t));
    b = zeros(1,length(t));
    c = zeros(1,length(t));
    d = zeros(1,length(t));

    a(1) = 0;
    b(1) = 0;
    c(1) = 0;
    d(1) = 0;

    for j=1:length(t)-1
        kaux1a = k1-k2*a(j)-k3*a(j)*b(j);
        kaux1b = k2*a(j)-k3*a(j)*b(j);
        kaux1c = k3*a(j)*b(j)-k4*c(j);
        kaux1d = k3*a(j)*b(j)-k5*d(j);
        auxa = a(j)+kaux1a*dt;
        auxb = b(j)+kaux1b*dt;
        auxc = c(j)+kaux1c*dt;
        auxd = d(j)+kaux1d*dt;
        kaux2a = k1-k2*auxa-k3*auxa*auxb;
        kaux2b = k2*auxa-k3*auxa*auxb;
        kaux2c = k3*auxa*auxb-k4*auxc;
        kaux2d = k3*auxa*auxb-k5*auxd;
        a(j+1) = a(j)+(kaux1a+kaux2a)*dt/2;
        b(j+1) = b(j)+(kaux1b+kaux2b)*dt/2;
        c(j+1) = c(j)+(kaux1c+kaux2c)*dt/2;
        d(j+1) = d(j)+(kaux1d+kaux2d)*dt/2;        
    end

    figure
    hold on
    plot(-100,-100,'-b','linewidth',2);
    plot(-100,-100,'-r','linewidth',2);
    plot(-100,-100,'-g','linewidth',2);
    plot(-100,-100,'-','Color',lightblueish,'linewidth',2);
    plot([0 Tmax],[0 0],'-k');
    plot(t,a,'-b','linewidth',2);
    plot(t,b,'-r','linewidth',2);
    plot(t,c,'-g','linewidth',2);
    plot(t,d,'-','Color',lightblueish,'linewidth',2);
    axis([0 Tmax -0.2 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [sec]');
    ylabel('Conc  [mM]');
    legend('a','b','c','d');
    %text(5,1.12,'k_1=k_4=k_5=0','fontsize',24)
    %text(5,1.12,'k_4=2  k_5=3','fontsize',24)
    
elseif MDL == 6
    
    % Parameters
    
    k1 = 9;
    k1m = 12;
    k2 = 2;
    
    k1 = 20;
    k1m = 3;
    k2 = 1;
    
    % Time definitions

    Tmax = 10;
    dt = 0.0001;
    t = 0:dt:Tmax;

    a = zeros(1,length(t));
    b = zeros(1,length(t));

    a(1) = 1;
    b(1) = 0;

    for j=1:length(t)-1
        kaux1a = -k1*a(j)+k1m*b(j);
        kaux1b = k1*a(j)-k1m*b(j)-k2*b(j);
        auxa = a(j)+kaux1a*dt;
        auxb = b(j)+kaux1b*dt;
        kaux2a = -k1*auxa+k1m*auxb;
        kaux2b = k1*auxa-k1m*auxb-k2*auxb;
        a(j+1) = a(j)+(kaux1a+kaux2a)*dt/2;
        b(j+1) = b(j)+(kaux1b+kaux2b)*dt/2;    
    end
    
    figure
    hold on
    plot(-100,-100,'-b','linewidth',2);
    plot(-100,-100,'-r','linewidth',2);
     plot(-100,-100,'-','Color',lightblueish,'linewidth',2);
    plot(t,a,'-b','linewidth',2);
    plot(t,b,'-r','linewidth',2);
     axis([-0.1 Tmax -0.2 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [sec]');
    ylabel('Conc  [mM]');
    legend('a','b');

    
end