% Solution of the Michaelis-Menten kinetic model
% Numerical method: modified Euler (Runge-Kutta, order 2)

clear all
close all

% ODEs:
% ds/dt = -s + (s+kpa-lda)c
% eps dc/dt = s - (s+kpa)c

KSE = 1;
        % 1: Dimensionless (2D)
        % 2: Dimensional (4D)

if KSE == 1        
        
    % Parameters
    
    kpa = 1;
    lda = 1;
    eps = 0.01;

    % Functions

    f = @(s,c) -s+(s+kpa-lda)*c;
    g = @(s,c) (s-(s+kpa)*c)/eps;


    Tmax = 10;
    dt = 0.005;
    t = 0:dt:Tmax;

    s = zeros(1,length(t));
    c = zeros(1,length(t));

    s(1) = 1;
    c(1) = 0;

    for j=1:length(t)-1
        k1s = f(s(j),c(j));
        k1c = g(s(j),c(j));
        as = s(j)+k1s*dt;
        ac = c(j)+k1c*dt;
        k2s = f(as,ac);   
        k2c = g(as,ac);
        s(j+1)=s(j)+(k1s+k2s)*dt/2;
        c(j+1)=c(j)+(k1c+k2c)*dt/2;    
    end
    
    ds = 0.01;
    svec = -1:ds:2;
    snlc = svec./(svec+kpa-lda);
    cnlc = svec./(svec+kpa);



    figure
    hold on
    plot(-100,-100,'-b','linewidth',2);
    plot(-100,-100,'-r','linewidth',2);
    plot([0 0],[-100 100],'--','Color',[.7 .7 .7],'linewidth',1.5);
    plot(t,s,'-b','linewidth',2);
    plot(t,c,'-r','linewidth',2);
    axis([-1 Tmax -0.1 1.1]);
    set(gca,'fontsize',20);
    xlabel('t');
    ylabel('');
    legend('S','C');
    
%     
%     
%     figure
%     hold on
%     plot(-100,-100,'-r','linewidth',2);
%     plot(-100,-100,'-g','linewidth',2);
%     plot(-100,-100,'-b','linewidth',2);
%     plot([0 0],[-10 10],'--','Color',[.7 .7 .7],'linewidth',1.5);
%     plot([-10 10],[0 0],'--','Color',[.7 .7 .7],'linewidth',1.5);
%     plot(svec,snlc,'-r','linewidth',2);
%     plot(svec,cnlc,'-g','linewidth',2);
%     plot(s,c,'-b','linewidth',2);
%     plot(s(1),c(1),'ob','linewidth',3);
%     axis([-0.5 2 -0.5 2]);
%     set(gca,'fontsize',20);
%     xlabel('S');
%     ylabel('C');
%     legend('s-nullcline','c-nullcline','trajectory');
    
MVN = 0;
if MVN == 1
    figure(101)
    k = 4;
    hFig = figure(101);
    set(hFig, 'Position', [20 100 1400 1000])
    for j=1:floor(length(t)/k)-1
        subplot(1,2,1)
        plot(-100,-100,'-b','linewidth',2);
        hold on
        plot(-100,-100,'-r','linewidth',2);
        plot([0 0],[-100 100],'--','Color',[.7 .7 .7],'linewidth',1.5);
        plot(t,s,'-b','linewidth',2);
        plot(t,c,'-r','linewidth',2);
        plot(t(k*j),s(k*j),'ob','linewidth',4);
        plot(t(k*j),c(k*j),'or','linewidth',4);
        axis([-1 Tmax -0.1 1.1]);
        set(gca,'fontsize',20);
        xlabel('t');
        ylabel('');
        legend('S','C');
        pause(0.0001);
        hold off  
        subplot(1,2,2)
        plot(-100,-100,'-r','linewidth',2);
         hold on
        plot(-100,-100,'-g','linewidth',2);
        plot(-100,-100,'-b','linewidth',2);
        plot([0 0],[-10 10],'--','Color',[.7 .7 .7],'linewidth',1.5);
        plot([-10 10],[0 0],'--','Color',[.7 .7 .7],'linewidth',1.5);
        plot(svec,snlc,'-r','linewidth',2);
        plot(svec,cnlc,'-g','linewidth',2);        
        plot(s(1),c(1),'ob','linewidth',3);
        plot(s(k*j),c(k*j),'ob','linewidth',3);
        plot(s(1:k*j),c(1:k*j),'--b','linewidth',1);
        axis([-0.5 2 -0.5 2]);
        set(gca,'fontsize',20);
        xlabel('S');
        ylabel('C');
        legend('s-nullcline','c-nullcline','trajectory');
        pause(0.0001);
        hold off  
    end
end
    
elseif KSE == 2
    
    % Parameters
    
    kp1 = 40;
    km1 = 1;
    k2 = 20;
    s0 = 5;
    c0 = 0;
    e0 = 0.5;    
    p0 = 0;
    
    % Functions
    
    fs = @(s,c,e)  -kp1*s*e+km1*c;
    fc = @(s,c,e)  -km1*c+kp1*s*e-k2*c;
    fe = @(s,c,e)  km1*c-kp1*s*e+k2*c;
    fp = @(c)      k2*c;
    
    Tmax = 5;
    dt = 0.01;
    t = 0:dt:Tmax;

    s = zeros(1,length(t));
    c = zeros(1,length(t));
    e = zeros(1,length(t));
    p = zeros(1,length(t));

    s(1) = s0;
    c(1) = c0;
    e(1) = e0;
    p(1) = p0;
    
    for j=1:length(t)-1
        k1s = fs(s(j),c(j),e(j));
        k1c = fc(s(j),c(j),e(j));
        k1e = fe(s(j),c(j),e(j));
        k1p = fp(c(j));        
        as = s(j)+k1s*dt;
        ac = c(j)+k1c*dt;
        ae = e(j)+k1e*dt;
        ap = p(j)+k1p*dt;
        k2s = fs(as,ac,ae);
        k2c = fc(as,ac,ae);
        k2e = fe(as,ac,ae);
        k2p = fp(ac);          
        s(j+1)=s(j)+(k1s+k2s)*dt/2;
        c(j+1)=c(j)+(k1c+k2c)*dt/2;
        e(j+1)=e(j)+(k1e+k2e)*dt/2;
        p(j+1)=p(j)+(k1p+k2p)*dt/2;
    end
    
    Ks = km1/kp1;
    Km = (km1+k2)/kp1;
    V = k2*c;
    Veq = k2*e0*s./(Ks+s);
    Vqss = k2*e0*s./(Km+s);
    
    figure
    hold on
    plot(-100,-100,'-b','linewidth',2);
    plot(-100,-100,'-r','linewidth',2);
    plot(-100,-100,'-g','linewidth',2);
    plot(-100,-100,'-c','linewidth',2);
    plot([0 0],[0 100],'--','Color',[.7 .7 .7],'linewidth',1.5);
    plot(t,s,'-b','linewidth',2);
    plot(t,c,'-r','linewidth',2);
    plot(t,e,'-g','linewidth',2);
    plot(t,p,'-c','linewidth',2);
    axis([-0.1 Tmax -0.1 6]);
    set(gca,'fontsize',20);
    xlabel('t');
    ylabel('');
    legend('S','C','E','P');
    
    
    figure
    hold on
    plot(-100,-100,'-b','linewidth',2);
    plot(-100,-100,'-r','linewidth',2);
    plot(-100,-100,'-g','linewidth',2);
    plot([0 0],[0 100],'--','Color',[.7 .7 .7],'linewidth',1.5);
    axis([0 Tmax -0.1 11]);
    plot(t,V,'-b','linewidth',2);
    plot(t,Veq,'-r','linewidth',2);
    plot(t,Vqss,'-g','linewidth',2);
    set(gca,'fontsize',20);
    xlabel('t');
    ylabel('');
    legend('V','V_{eq}','V_{qss}');
    
end
