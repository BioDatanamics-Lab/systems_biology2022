% Covalent modification cycle & cascades
% 2022-04-26

clearvars;
close all;

% S + E1 <--> C1 ---> P + E1
% P + E2 <--> C2 ---> S + E2
% ReactI0n constants: ksc1, kcs1, kcp1, kpc2, kcp2, kcs2
% 
% NotatI0n:
%           kxy: S-->C reactI0n (irreversible network component)
% 

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

% Functions 

heaviside=@(t) 0.5*(t == 0)+(t > 0);


MDL = 1;
        % 1: Covalent modificatI0n cycle
        
if MDL == 1

    % Parameters
    
    ksc1 = 100;
    kcs1 = 10;
    kcp1 = 2;
    kpc2 = 10;
    kcp2 = 10;
    kcs2 = 2;
    
    S0 = 0;
    C10 = 0;
    C20 = 0;
    E10 = 0.25;
    E20 = 0.25;
    P0 = 0;

    E1t = E10;
    E2t = E20;
    
    Tmax = 100;
    dt = 0.01;
    t = 0:dt:Tmax;
    
    ton = 20;
    toff = 25;
    Ain1 = 0.1;
    Fsqw = Ain1*heaviside(t-ton).*heaviside(toff-t);
    
    S = zeros(1,length(t));
    E1 = zeros(1,length(t));
    E2 = zeros(1,length(t));
    P = zeros(1,length(t));
    
    S(1) = S0;
    E1(1) = E10;
    E2(1) = E20;
    P(1) = P0;
    
    for j=1:length(t)-1
        k1s = -ksc1*S(j)*E1(j)+kcs1*(E1t-E1(j))+kcs2*(E2t-E2(j));
        k1s = k1s+Fsqw(j);
        k1e1 = -ksc1*S(j)*E1(j)+(kcs1+kcp1)*(E1t-E1(j));
        k1e2 = -kpc2*P(j)*E2(j)+(kcp2+kcs2)*(E2t-E2(j));
        k1p = kcp1*(E1t-E1(j))-kpc2*P(j)*E2(j)+kcp2*(E2t-E2(j));
        as = S(j)+k1s*dt;
        ae1 = E1(j)+k1e1*dt;
        ae2 = E2(j)+k1e2*dt;
        ap = P(j)+k1p*dt;     
        k2s = -ksc1*S(j)*ae1+kcs1*(E1t-ae1)+kcs2*(E2t-ae2);
        k2s = k2s+Fsqw(j);
        k2e1 = -ksc1*S(j)*ae1+(kcs1+kcp1)*(E1t-ae1);
        k2e2 = -kpc2*P(j)*ae2+(kcp2+kcs2)*(E2t-ae2);
        k2p = kcp1*(E1t-ae1)-kpc2*P(j)*ae2+kcp2*(E2t-ae2);
        S(j+1) = S(j)+(k1s+k2s)*dt/2;
        E1(j+1) = E1(j)+(k1e1+k2e1)*dt/2;
        E2(j+1) = E2(j)+(k1e2+k2e2)*dt/2;
        P(j+1) = P(j)+(k1p+k2p)*dt/2;
    end
    
    C1 = E1t-E1;
    C2 = E2t-E2;
   
    
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
    %plot(-100,-100,'-r','linewidth',2);
    %plot(-100,-100,'-g','linewidth',2);
    plot(-100,-100,'-','Color',lightblueish,'linewidth',2);
    %plot(-100,-100,'-','Color',lightcoral,'linewidth',2);
    %plot(-100,-100,'-','Color',mediumacquamarine,'linewidth',2);
    plot([0 Tmax],[0 0],'-k');
    plot([0 0],[-10 10],'-k');
    plot(t,S,'-b','linewidth',2);
    %plot(t,C1,'-r','linewidth',2);
    %plot(t,E1,'-g','linewidth',2);
    plot(t,P,'-','Color',lightblueish,'linewidth',2);
    %plot(t,C2,'-','Color',lightcoral,'linewidth',2);
    %plot(t,E2,'-','Color',mediumacquamarine,'linewidth',2);
    axis([-0.05 Tmax -0.2 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [au]');
    ylabel('Conc  [au]');
    %legend('S','C_1','P','C_2');
    legend('S','P');

end