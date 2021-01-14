% Feedback controllers for MNB Dist Rej
% ============================================================
%   Copyright (c) 2008-, Xu Chen, chx@uw.edu
%   Author(s): Xu Chen
% ============================================================
folder = fileparts(which(mfilename));
addpath(genpath(folder));

bodeopt = bodeoptions;
bodeopt.FreqUnits = 'Hz';
bodeopt.XLim = [1 0.5/Ts];
bode_opt = bodeopt;

Notch_w         = [2*pi*4100 2*pi*8200 2*pi*12300 2*pi*5000];
Notch_Zeta_den  = [0.394, 0.22,  0.47, 0.002];
Notch_Zeta_num  = [0.02,  0.02,  0.02, 0.001];

NF_4100     = tf([1 2*Notch_Zeta_num(1)*Notch_w(1) Notch_w(1)^2],...
    [1 2*Notch_Zeta_den(1)*Notch_w(1) Notch_w(1)^2]);
sys_NF_4100 = c2d(NF_4100,Ts,'matched');

NF_8200     = tf([1 2*Notch_Zeta_num(2)*Notch_w(2) Notch_w(2)^2],...
    [1 2*Notch_Zeta_den(2)*Notch_w(2) Notch_w(2)^2]);
sys_NF_8200 = c2d(NF_8200,Ts,'matched');

NF_12300    = tf([1 2*Notch_Zeta_num(3)*Notch_w(3) Notch_w(3)^2],...
    [1 2*Notch_Zeta_den(3)*Notch_w(3) Notch_w(3)^2]);
sys_NF_12300 = c2d(NF_12300,Ts,'matched');

NF_5000     = tf([1 2*Notch_Zeta_num(4)*Notch_w(4) Notch_w(4)^2],...
    [1 2*Notch_Zeta_den(4)*Notch_w(4) Notch_w(4)^2]);
sys_NF_5000 = c2d(NF_5000,Ts,'matched');

% Series connection of notch filters
sys_NF      = sys_NF_4100 * sys_NF_8200 * sys_NF_12300 * sys_NF_5000;

% PID controller
Freq_zc  = 1000;                 % Hz
coef_PID = [2,0.05,70];          % P,I,D gain
sys_C    = tf(coef_PID(1)*[1,-1,0]+...
    coef_PID(2)*[0,1,0]+...
    coef_PID(3)*[1,-2,1],...
    [1, -1, 0], Ts);

% Gain adjustment
g_adjust = 1.665e3;
sys_C    = sys_C/g_adjust;

robust_sw = 2; % default
% Switch for robustness test

[NmaxTs,m] = size(SensorNoise.Data);
[NmaxTc,m] = size(ForceDist.Data);
sys_C1 = sys_C;
sys_C2 = sys_NF;
%%
if 0
    NF2x_4100     = tf([1 2*Notch_Zeta_num(1)*Notch_w(1) Notch_w(1)^2],...
        [1 2*Notch_Zeta_den(1)*Notch_w(1) Notch_w(1)^2]);
    sys_NF2x_4100 = c2d(NF_4100,Ts/2,'matched');
    
    NF2x_8200     = tf([1 2*Notch_Zeta_num(2)*Notch_w(2) Notch_w(2)^2],...
        [1 2*Notch_Zeta_den(2)*Notch_w(2) Notch_w(2)^2]);
    sys_NF2x_8200 = c2d(NF_8200,Ts/2,'matched');
    
    NF2x_12300    = tf([1 2*Notch_Zeta_num(3)*Notch_w(3) Notch_w(3)^2],...
        [1 2*Notch_Zeta_den(3)*Notch_w(3) Notch_w(3)^2]);
    sys_NF2x_12300 = c2d(NF_12300,Ts/2,'matched');
    
    NF2x_5000     = tf([1 2*Notch_Zeta_num(4)*Notch_w(4) Notch_w(4)^2],...
        [1 2*Notch_Zeta_den(4)*Notch_w(4) Notch_w(4)^2]);
    sys_NF2x_5000 = c2d(NF_5000,Ts/2,'matched');
    
    % Series connection of notch filters
    sys_NF2x      = sys_NF2x_4100 * sys_NF2x_8200 * sys_NF2x_12300 * sys_NF2x_5000;
    
    NF4x_4100     = tf([1 2*Notch_Zeta_num(1)*Notch_w(1) Notch_w(1)^2],...
        [1 2*Notch_Zeta_den(1)*Notch_w(1) Notch_w(1)^2]);
    sys_NF4x_4100 = c2d(NF_4100,Ts/4,'matched');
    
    NF4x_8200     = tf([1 2*Notch_Zeta_num(2)*Notch_w(2) Notch_w(2)^2],...
        [1 2*Notch_Zeta_den(2)*Notch_w(2) Notch_w(2)^2]);
    sys_NF4x_8200 = c2d(NF_8200,Ts/4,'matched');
    
    NF4x_12300    = tf([1 2*Notch_Zeta_num(3)*Notch_w(3) Notch_w(3)^2],...
        [1 2*Notch_Zeta_den(3)*Notch_w(3) Notch_w(3)^2]);
    sys_NF4x_12300 = c2d(NF_12300,Ts/4,'matched');
    
    NF4x_5000     = tf([1 2*Notch_Zeta_num(4)*Notch_w(4) Notch_w(4)^2],...
        [1 2*Notch_Zeta_den(4)*Notch_w(4) Notch_w(4)^2]);
    sys_NF4x_5000 = c2d(NF_5000,Ts/4,'matched');
    
    % Series connection of notch filters
    sys_NF4x      = sys_NF4x_4100 * sys_NF4x_8200 * sys_NF4x_12300 * sys_NF4x_5000;
    
    figure,
    bodeplot(sys_C2,sys_NF2x,'r--',sys_NF4x,'m:',bode_opt)
    %     hold on
    legend('1x NF','2x NF','4x NF')
    h = findobj(gcf,'type','line');
    set(h,'linewidth',1.5);
end
%%
% save benchmark_feedback_controller sys_C1 sys_C2
