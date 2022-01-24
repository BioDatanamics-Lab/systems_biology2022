% 2022-01-23
% Original code: 2020-10-20


% SIR Model: basic

% S: Suscetible
% I: Infectious
% R: Removed

clearvars;
close all;

MDL = 3;
%       1: SIR model 
%       2: SIR model + negative feedback
%       3: SIR model + negative feedback + positive feedback (quadratic)
%
% Note: 3 reduces to 1 (lda=gma=0) and 2 (gma=0)

if MDL == 1
    % Parameters

    a = 0.5;
    c = 0.2;
    N = 100;


    % Time definitions

    Tmax = 100;
    dt = 0.01;
    t = 0:dt:Tmax;

    % Numerical Simulations

    S = zeros(1,length(t));
    I = zeros(1,length(t));
    R = zeros(1,length(t));

    S(1) = 99;
    I(1) = N-S(1);
    R(1) = 0;

    for j=1:length(t)-1
        k1s = -a*S(j)*I(j)/N;
        k1i = a*S(j)*I(j)/N-c*I(j);
        k1r = c*I(j);
        as = S(j)+k1s*dt;
        ai = I(j)+k1i*dt;
        ar = R(j)+k1r*dt;
        k2s = -a*as*ai/N;
        k2i = a*as*ai/N-c*ai;
        k2r = c*ai;
        S(j+1) = S(j)+(k1s+k2s)*dt/2;
        I(j+1) = I(j)+(k1i+k2i)*dt/2;
        R(j+1) = R(j)+(k1r+k2r)*dt/2;
    end

    figure
    hold on
    plot(t,S,'--b','linewidth',2);
    plot(t,I,'-r','linewidth',2);
    plot(t,R,'--g','linewidth',2);
    set(gca,'fontsize',24);
    xlabel('t   [au]');
    legend('S','I','R');

elseif MDL == 2
    
    % Parameters

    a = 0.5;
    c = 0.2;
    lda = 0.01;
    N = 100;


    % Time definitions

    Tmax = 400;
    dt = 0.01;
    t = 0:dt:Tmax;

    % Numerical Simulations

    S = zeros(1,length(t));
    I = zeros(1,length(t));
    R = zeros(1,length(t));

    S(1) = 99;
    I(1) = N-S(1);
    R(1) = 0;

    for j=1:length(t)-1
        k1s = -a*S(j)*I(j)/N+lda*R(j);
        k1i = a*S(j)*I(j)/N-c*I(j);
        k1r = c*I(j)-lda*R(j);
        as = S(j)+k1s*dt;
        ai = I(j)+k1i*dt;
        ar = R(j)+k1r*dt;
        k2s = -a*as*ai/N+lda*ar;
        k2i = a*as*ai/N-c*ai;
        k2r = c*ai-lda*ar;
        S(j+1) = S(j)+(k1s+k2s)*dt/2;
        I(j+1) = I(j)+(k1i+k2i)*dt/2;
        R(j+1) = R(j)+(k1r+k2r)*dt/2;
    end

    figure
    hold on
    plot(t,S,'--b','linewidth',2);
    plot(t,I,'-r','linewidth',2);
    plot(t,R,'--g','linewidth',2);
    set(gca,'fontsize',24);
    xlabel('t   [au]');
    legend('S','I','R');

elseif MDL == 3

    % Parameters

    a = 0.5;
    c = 0.2;
    lda = 0.01;
    gma = 0.1;
    N = 100;

    % Time definitions

    Tmax = 1000;
    dt = 0.01;
    t = 0:dt:Tmax;

    % Numerical Simulations

    S = zeros(1,length(t));
    I = zeros(1,length(t));
    R = zeros(1,length(t));

    S(1) = 99;
    I(1) = N-S(1);
    R(1) = 0;

    for j=1:length(t)-1
        k1s = -a*(1+gma*I(j))*S(j)*I(j)/N+lda*R(j);
        k1i = a*(1+gma*I(j))*S(j)*I(j)/N-c*I(j);
        k1r = c*I(j)-lda*R(j);
        as = S(j)+k1s*dt;
        ai = I(j)+k1i*dt;
        ar = R(j)+k1r*dt;
        k2s = -a*(1+gma*ai)*as*ai/N+lda*ar;
        k2i = a*(1+gma*ai)*as*ai/N-c*ai;
        k2r = c*ai-lda*ar;
        S(j+1) = S(j)+(k1s+k2s)*dt/2;
        I(j+1) = I(j)+(k1i+k2i)*dt/2;
        R(j+1) = R(j)+(k1r+k2r)*dt/2;
    end

    figure
    hold on
    plot(t,S,'--b','linewidth',2);
    plot(t,I,'-r','linewidth',2);
    plot(t,R,'--g','linewidth',2);
    set(gca,'fontsize',24);
    xlabel('t   [au]');
    legend('S','I','R');

end




