% 2022-02-11

% Enzymatic reactI0ns and covalent modificatI0n cycles

% Enzymatic reactI0n
%
% S + E <--> C ---> P + E
% ReactI0n constants: ksc, kcs & kcp
%
% Covalent modificatI0n cycle
%
% S + E1 <--> C1 ---> P + E1
% P + E2 <--> C2 ---> S + E2
% ReactI0n constants: ksc1, kcs1, kcp1, kpc2, kcp2, kcs2
% 
% NotatI0n:
%           kxy: S-->C reactI0n (irreversible network component)
% 

clearvars;
close all;

lightblueish = [.4 .6 .9];
lightcoral = [0.94 0.5 0.5];
lightsalmon = [0.9 0.5 0.4];
lightgray = [.7 .7 .7];
darkgray = [.3 .3 .3];
darkgray2 = [.1 .1 .1];
mediumacquamarine = [0.4 0.8 0.6];

% [t] = au
% [conc] = au
% ksc = conc^-1 t^-1
% kcs = t^-1
% kcp = t^-1

MDL = 2;
        % 1: Enzymatic reaction
        % 2: Covalent modificatI0n cycle
        % 3: Competitive inhibitI0n
        % 4: Allosteric inhibitI0n
        % 5: Cooperativity
        
        
if MDL == 1
    
    % Parameters
    
%     ksc = 30;
%     kcs = 1;
%     kcp = 10;
%     S0 = 5;
%     E0 = 1; 
    
    ksc = 100;
    kcs = 10;
    kcp = 2;
    
    S0 = 1;
    C0 = 0;
    E0 = 0.2;    
    P0 = 0;

    Et = E0;
    
    Tmax = 10;
    dt = 0.001;
    t = 0:dt:Tmax;
    
    S = zeros(1,length(t));
    C = zeros(1,length(t));
    
    S(1) = S0;
    C(1) = C0;
   
    for j=1:length(t)-1
        k1s = -ksc*S(j)*(Et-C(j))+kcs*C(j);
        k1c = ksc*S(j)*(Et-C(j))-(kcs+kcp)*C(j);
        as = S(j)+k1s*dt;
        ac = C(j)+k1c*dt;
        k2s = -ksc*as*(Et-ac)+kcs*ac;
        k2c = ksc*as*(Et-ac)-(kcs+kcp)*ac;
        S(j+1) = S(j)+(k1s+k2s)*dt/2;
        C(j+1) = C(j)+(k1c+k2c)*dt/2;        
    end
    
    E = Et-C;
    P = S0+P0-S-C;
    
    figure
    hold on
    plot(-100,-100,'-b','linewidth',2);
    plot(-100,-100,'-r','linewidth',2);
    plot(-100,-100,'-g','linewidth',2);
    plot(-100,-100,'-','Color',lightblueish,'linewidth',2);
    plot([0 Tmax],[0 0],'-k');
    plot([0 0],[-10 10],'-k');
    plot(t,S,'-b','linewidth',2);
    plot(t,C,'-r','linewidth',2);
    plot(t,E,'-g','linewidth',2);
    plot(t,P,'-','Color',lightblueish,'linewidth',2);
    axis([-0.05 Tmax -0.2 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [au]');
    ylabel('Conc  [au]');
    legend('S','C','E','P');
    
elseif MDL == 2
    
    % Parameters
    
    ksc1 = 100;
    kcs1 = 10;
    kcp1 = 2;
    kpc2 = 100;
    kcp2 = 10;
    kcs2 = 2;
    
    S0 = 1;
    C10 = 0;
    C20 = 0;
    E10 = 0.25;
    E20 = 0.25;
    P0 = 0;

    E1t = E10;
    E2t = E20;
    
    Tmax = 10;
    dt = 0.01;
    t = 0:dt:Tmax;
    
    C1 = zeros(1,length(t));
    C2 = zeros(1,length(t));
    P = zeros(1,length(t));
    
    C1(1) = C10;
    C2(1) = C20;
    P(1) = P0;
    
    for j=1:length(t)-1
        k1c1 = ksc1*(S0+P0-P(j)-C1(j)-C2(j))*(E1t-C1(j))-(kcs1+kcp1)*C1(j);
        k1c2 = kpc2*P(j)*(E2t-C2(j))-(kcp2+kcs2)*C2(j);
        k1p = kcp1*C1(j)-kpc2*P(j)*(E2t-C2(j))+kcp2*C2(j);
        ac1 = C1(j)+k1c1*dt;
        ac2 = C2(j)+k1c2*dt;
        ap = P(j)+k1p*dt;
        k2c1 = ksc1*(S0+P0-ap-ac1-ac2)*(E1t-ac1)-(kcs1+kcp1)*ac1;
        k2c2 = kpc2*ap*(E2t-ac2)-(kcp2+kcs2)*ac2;
        k2p = kcp1*ac1-kpc2*ap*(E2t-ac2)+kcp2*ac2;
        C1(j+1) = C1(j)+(k1c1+k2c1)*dt/2;
        C2(j+1) = C2(j)+(k1c2+k2c2)*dt/2;
        P(j+1) = P(j)+(k1p+k2p)*dt/2;
    end
    
    E1 = E1t-C1;
    E2 = E2t-C2;
    S = S0+P0-P-C1-C2;
    
    figure
    hold on
    plot(-100,-100,'-b','linewidth',2);
    plot(-100,-100,'-r','linewidth',2);
    plot(-100,-100,'-g','linewidth',2);
    plot(-100,-100,'-','Color',lightblueish,'linewidth',2);
    plot(-100,-100,'-','Color',lightcoral,'linewidth',2);
    plot(-100,-100,'-','Color',mediumacquamarine,'linewidth',2);
    plot([0 Tmax],[0 0],'-k');
    plot([0 0],[-10 10],'-k');
    plot(t,S,'-b','linewidth',2);
    plot(t,C1,'-r','linewidth',2);
    plot(t,E1,'-g','linewidth',2);
    plot(t,P,'-','Color',lightblueish,'linewidth',2);
    plot(t,C2,'-','Color',lightcoral,'linewidth',2);
    plot(t,E2,'-','Color',mediumacquamarine,'linewidth',2);
    axis([-0.05 Tmax -0.2 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [au]');
    ylabel('Conc  [au]');
    legend('S','C_1','E_1','P','C_2','E_2');
    
    figure
    hold on
    plot(-100,-100,'-b','linewidth',2);
    plot(-100,-100,'-r','linewidth',2);
    %plot(-100,-100,'-g','linewidth',2);
    plot(-100,-100,'-','Color',lightblueish,'linewidth',2);
    plot(-100,-100,'-','Color',lightcoral,'linewidth',2);
    %plot(-100,-100,'-','Color',mediumacquamarine,'linewidth',2);
    plot([0 Tmax],[0 0],'-k');
    plot([0 0],[-10 10],'-k');
    plot(t,S,'-b','linewidth',2);
    plot(t,C1,'-r','linewidth',2);
    %plot(t,E1,'-g','linewidth',2);
    plot(t,P,'-','Color',lightblueish,'linewidth',2);
    plot(t,C2,'-','Color',lightcoral,'linewidth',2);
    %plot(t,E2,'-','Color',mediumacquamarine,'linewidth',2);
    axis([-0.05 Tmax -0.2 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [au]');
    ylabel('Conc  [au]');
    legend('S','C_1','P','C_2');
    
elseif MDL == 3
    
    ksc = 100;
    kcs = 10;
    kcp = 2;
    kid = 100;
    kdi = 0;
    
    S0 = 1;
    C0 = 0;
    E0 = 0.2;    
    P0 = 0;
    I0 = 1;
    D0 = 0;

    Et = E0;
    
    Tmax = 6;
    dt = 0.001;
    t = 0:dt:Tmax;
    
    S = zeros(1,length(t));
    C = zeros(1,length(t));
    D = zeros(1,length(t));
    
    S(1) = S0;
    C(1) = C0;
    D(1) = D0;
   
    for j=1:length(t)-1
        k1s = -ksc*S(j)*(Et-C(j)-D(j))+kcs*C(j);
        k1c = ksc*S(j)*(Et-C(j)-D(j))-(kcs+kcp)*C(j);
        k1d = kid*(I0-D(j))*(Et-C(j)-D(j))-kdi*D(j);
        as = S(j)+k1s*dt;
        ac = C(j)+k1c*dt;
        ad = D(j)+k1d*dt;
        k2s = -ksc*as*(Et-ac)+kcs*ac;
        k2c = ksc*as*(Et-ac-ad)-(kcs+kcp)*ac;
        k2d = kid*(I0-ad)*(Et-ac-ad)-kdi*ad;
        S(j+1) = S(j)+(k1s+k2s)*dt/2;
        C(j+1) = C(j)+(k1c+k2c)*dt/2;       
        D(j+1) = D(j)+(k1d+k2d)*dt/2;
    end
    
    
    E = Et-C-D;
    P = S0+P0-S-C;
    I = I0-D;
   
    
    figure
    hold on
    plot(-100,-100,'-b','linewidth',2);
    plot(-100,-100,'-r','linewidth',2);
    plot(-100,-100,'-g','linewidth',2);
    plot(-100,-100,'-','Color',lightblueish,'linewidth',2);
    plot([0 Tmax],[0 0],'-k');
    plot([0 0],[-10 10],'-k');
    plot(t,S,'-b','linewidth',2);
    plot(t,C,'-r','linewidth',2);
    plot(t,E,'-g','linewidth',2);
    plot(t,P,'-','Color',lightblueish,'linewidth',2);
    axis([-0.05 Tmax -0.2 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [au]');
    ylabel('Conc  [au]');
    legend('S','C','E','P');
    
    figure
    hold on
    plot(-100,-100,'--b','linewidth',2);
    plot(-100,-100,'--r','linewidth',2);
    plot(-100,-100,'--g','linewidth',2);
    plot(-100,-100,'--','Color',lightblueish,'linewidth',2);
    plot([0 Tmax],[0 0],'-k');
    plot([0 0],[-10 10],'-k');
    plot(t,I,'--b','linewidth',2);
    plot(t,D,'--r','linewidth',2);
    axis([-0.05 Tmax -0.2 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [au]');
    ylabel('Conc  [au]');
    legend('I','D');
    
elseif MDL == 4
    
    ksc = 100;
    kcs = 100;
    kcp = 2;
    ksf = 0;
    kfs = 0;
    kif = 0;
    kfi = 0;
    kid = 0.1;
    kdi = 0;
    
    S0 = 1;
    C0 = 0;
    E0 = 0.2;    
    P0 = 0;
    I0 = 1;
    D0 = 0;
    F0 = 0;
    
    Et = E0;
    
    Tmax = 6;
    dt = 0.001;
    t = 0:dt:Tmax;
    
    S = zeros(1,length(t));
    C = zeros(1,length(t));
    D = zeros(1,length(t));
    F = zeros(1,length(t));
    
    S(1) = S0;
    C(1) = C0;
    D(1) = D0;
    F(1) = F0;
   
    for j=1:length(t)-1
        k1s = -ksc*S(j)*(Et-C(j)-D(j)-F(j))+kcs*C(j)-ksf*S(j)*D(j)+kfs*F(j);
        k1c = ksc*S(j)*(Et-C(j)-D(j)-F(j))-(kcs+kcp)*C(j)-kif*(I0-D(j)-F(j))+kfi*F(j);
        k1d = kid*(I0-D(j)-F(j))*(Et-C(j)-D(j)-F(j))-kdi*D(j)-ksf*S(j)*D(j)+kfs*F(j);
        k1f = kif*(I0-D(j)*C(j))-kfi*F(j)+ksf*S(j)*D(j)-kfs*F(j);
        as = S(j)+k1s*dt;
        ac = C(j)+k1c*dt;
        ad = D(j)+k1d*dt;
        af = F(j)+k1f*dt;
        k2s = -ksc*as*(Et-ac)+kcs*ac-ksf*as*ad+kfs*af;
        k2c = ksc*as*(Et-ac-ad-af)-(kcs+kcp)*ac-kif*(I0-ad-af)*ac+kfi*af;
        k2d = kid*(I0-ad-af)*(Et-ac-ad-af)-kdi*ad-ksf*as*ad+kfs*af;
        k2f = kif*(I0-ad-af)*ac-kfi*af+ksf*af+ksf*as*ad-kfs*af;
        S(j+1) = S(j)+(k1s+k2s)*dt/2;
        C(j+1) = C(j)+(k1c+k2c)*dt/2;       
        D(j+1) = D(j)+(k1d+k2d)*dt/2;
        F(j+1) = F(j)+(k1f+k2f)*dt/2;
    end
    
    E = Et-C-D-F;
    P = S0+P0-S-C-F;
    I = I0-D-F;
    
    figure
    hold on
    plot(-100,-100,'-b','linewidth',2);
    plot(-100,-100,'-r','linewidth',2);
    plot(-100,-100,'-g','linewidth',2);
    plot(-100,-100,'-','Color',lightblueish,'linewidth',2);
    plot([0 Tmax],[0 0],'-k');
    plot([0 0],[-10 10],'-k');
    plot(t,S,'-b','linewidth',2);
    plot(t,C,'-r','linewidth',2);
    plot(t,E,'-g','linewidth',2);
    plot(t,P,'-','Color',lightblueish,'linewidth',2);
    axis([-0.05 Tmax -0.2 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [au]');
    ylabel('Conc  [au]');
    legend('S','C','E','P');
    
    figure
    hold on
    plot(-100,-100,'--b','linewidth',2);
    plot(-100,-100,'--r','linewidth',2);
    plot(-100,-100,'--g','linewidth',2);
    plot(-100,-100,'--','Color',lightblueish,'linewidth',2);
    plot([0 Tmax],[0 0],'-k');
    plot([0 0],[-10 10],'-k');
    plot(t,I,'--b','linewidth',2);
    plot(t,D,'--r','linewidth',2);
    plot(t,F,'--g','linewidth',2);
    axis([-0.05 Tmax -0.2 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [au]');
    ylabel('Conc  [au]');
    legend('I','D','F');
        
elseif MDL == 5
    
   k1 = 1;
   km1 = 1;
   k2 = 1;
   km2 = 1;
   k3 = 1;
   km3 = 1;
   k4 = 1;
   km4 = 1;
   
   P0 = 1;
   X0 = 1;
   
   
    Tmax = 10;
    dt = 0.001;
    t = 0:dt:Tmax;
    
    P = zeros(1,length(t));
    P1 = zeros(1,length(t));
    P2 = zeros(1,length(t));
    P3 = zeros(1,length(t));
    P4 = zeros(1,length(t));
    X = zeros(1,length(t));
    
    P(1) = P0;
    X(1) = X0;
    
    for j=1:length(t)-1
        k1p = -4*k1*P(j)*X(j)+km1*P1(j);
        k1p1 = 4*k1*P(j)*X(j)-km1*P1(j)-3*k2*P1(j)*X(j)+2*km2*P2(j);
        k1p2 = 3*k2*P1(j)*X(j)-2*km2*P2(j)-2*k3*P2(j)*X(j)+3*km3*P3(j);
        k1p3 = 2*k3*P2(j)*X(j)-3*km3*P3(j)-k4*P3(j)*X(j)+4*km4*(P0-P(j)-P1(j)-P2(j)-P3(j));
        k1x = -(4*k1*P(j)+3*k2*P1(j)+2*k3*P2(j)+k4*P3(j))*X(j)+km1*P1(j)+2*km2*P2(j)+3*km3*P3(j)+4*km4*(P0-P(j)-P1(j)-P2(j)-P3(j));
        ap = P(j)+k1p*dt;
        ap1 = P1(j)+k1p1*dt;
        ap2 = P2(j)+k1p2*dt;
        ap3 = P3(j)+k1p3*dt;
        ax = X(j)+k1x*dt;
        k2p = -4*k1*ap*ax+km1*ap1;
        k2p1 = 4*k1*ap*ax-km1*ap1-3*k2*ap1*ax+2*km2*ap2;
        k2p2 = 3*k2*ap1*ax-2*km2*ap2-2*k3*ap2*ax+3*km3*ap3;
        k2p3 = 2*k3*ap2*ax-3*km3*ap3-k4*ap3*ax+4*km4*(P0-ap-ap1-ap2-ap3);
        k2x = -(4*k1*ap+3*k2*ap1+2*k3*ap2+k4*ap3)*ax+km1*ap1+2*km2*ap2+3*km3*ap3+4*km4*(P0-ap-ap1-ap2-ap3);
        P(j+1) = P(j)+(k1p+k2p)*dt/2;
        P1(j+1) = P1(j)+(k1p1+k2p1)*dt/2;
        P2(j+1) = P2(j)+(k1p2+k2p2)*dt/2;
        P3(j+1) = P3(j)+(k1p3+k2p3)*dt/2;
        X(j+1) = X(j)+(k1x+k2x)*dt/2;
    end
    
    figure
    hold on
    plot(-100,-100,'-b','linewidth',2);
    plot(-100,-100,'-r','linewidth',2);
    plot(-100,-100,'-g','linewidth',2);
    plot(-100,-100,'-','Color',lightblueish,'linewidth',2);
    plot(-100,-100,'-','Color',lightcoral,'linewidth',2);
    plot(-100,-100,'-','Color',mediumacquamarine,'linewidth',2);
    plot([0 Tmax],[0 0],'-k');
    plot([0 0],[-10 10],'-k');
    plot(t,P,'-b','linewidth',2);
    plot(t,P1,'-r','linewidth',2);
    plot(t,P2,'-g','linewidth',2);
    plot(t,P3,'-','Color',lightblueish,'linewidth',2);
    plot(t,P4,'-','Color',lightcoral,'linewidth',2);
    plot(t,X,'-','Color',mediumacquamarine,'linewidth',2);
    axis([-0.05 Tmax -0.2 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [au]');
    ylabel('Conc  [au]');
    legend('P','P_1','P_2','P_3','P_4','X');
    
    
   
end
        
        
      