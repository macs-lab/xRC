function Pn_ZPETinv = NominalPselect(Pn,sysd_P,Ts,method,plotflag)
% Nominal plant select
% function Pn_ZPETinv = NominalPselect(Pn,sysd_P,Ts,method,plotflag)
% 
% G_nominal inverse for DOB with delay shifted to the left of DOB
% [num_Pn, den_Pn] = ss2tf(Pn.a, Pn.b, Pn.c, Pn.d);
% inputs: 
%       Pn:     nominal transfer function in CT
%       sysd_P: full order transfer function in DT
% ============================================================
%   Copyright (c) 2008-, Xu Chen, xuchen@cal.berkeley.edu
%   Author(s): Xu Chen 
% ============================================================


% Pn_tf = tf(num_Pn, den_Pn);
% figure;bode(Pn,Pn_tf);legend('Pn','Pn tf');title('tf loses the delay')
% 1.137e-013 s + 3.745e009
% -------------------------
% s^2 + 565.5 s + 3.198e005

Pn_d = c2d(Pn,Ts,'zoh');
[num_Pn_d, den_Pn_d] = ss2tf(Pn_d.a, Pn_d.b, Pn_d.c, Pn_d.d);
disp 'P:'
Pn_d_tf = tf(num_Pn_d, den_Pn_d,Ts)
% figure;bode(Pn,Pn_d_tf);legend('Pn','Pn d tf')
% figure;bode(Pn_d_tf,Pn_tf);legend('Pn d tf','Pn tf')
% 1.448 z^2 + 3.685 z + 0.1836
% ----------------------------
%  z^3 - 1.978 z^2 + 0.9788 z
%
% Sampling time: 3.7879e-005

if strcmp(method,'zpet')
    Pn_ZPETinv = ZPETinv(Pn_d,Ts);%ZPET inverse of Pn_d
    % may have some problems 2011-09-18
elseif strcmp(method,'sInv')
    Pn_ZPETinv = stable_inv(Pn_d,Ts,'noplot',0.8,1,sysd_P);%(replace the unstable zero with a stable one)
%     function Pinv =
%     stable_inv(P,Ts,plotflag,z_replace,unstable_gard,P_full)
else % default selection
%     Pn_ZPETinv = stable_inv(Pn_d,Ts,'noplot');%(replace the unstable zero
%     with a stable one)
    Pn_ZPETinv = stable_inv(Pn_d,Ts,'noplot',0.8,1,sysd_P);%(replace the unstable zero with a stable one)
end
disp 'P_inv'
Pn_ZPETinv.tf

if strcmp(plotflag,'plot')
    % [LP,LP_c] = Qfilter(0.00008,Ts,1);%0.0005 0.00008
    % G_n_inv0 = z^-1/Pn_d_tf*LP;
    % G_n_inv0 = z^-1/Pn_d_tf;
    
    z = tf('z',Ts);
    
    figure;
    h = bodeplot(z^Pn_ZPETinv.delaySteps*Pn_ZPETinv.tf);
    setoptions(h, 'FreqUnits', 'Hz', 'Grid','on');
    xlim([1,0.5/Ts])
    hold on;
    bodeplot(z^Pn_ZPETinv.delaySteps*Pn_ZPETinv.tf * Pn_d_tf,'r--');
    bodeplot(z^Pn_ZPETinv.delaySteps*Pn_ZPETinv.tf * sysd_P,'k-');
%     ylim([-360,360])
    legend('Nominal plant stable inverse', 'Nominal P * its Inverse' , 'Inverse * full order P','Location','Best'); 
    
    
    figure,
    h = bodeplot(sysd_P,'r--');
    setoptions(h, 'FreqUnits', 'Hz', 'Grid','on','PhaseWrapping','On');
    h = findobj(gcf,'type','line');
    set(h,'linewidth',2);
    xlim([10,0.5/Ts])
    hold on;
    bodeplot(z^-Pn_ZPETinv.delaySteps/Pn_ZPETinv.tf);
    h = findobj(gcf,'type','line');
    set(h,'linewidth',1.5);
    legend('Full-order plant', 'Nominal model','Location','Best'); 
    
%     figure;pzmap(Pn_ZPETinv.tf);title 'ZPET Inverse of nominal plant';
%     figure;pzmap(Pn_d_tf);title 'DT Nominal plant';
%     figure;pzmap(Pn);title 'CT Pn';
%     figure;pzmap(1/Pn_tf);title '1/Pn';
end
