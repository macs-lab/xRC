function Pinv = stable_inv(P,Ts,plotflag,z_replace,unstable_gard,P_full)
% function Pinv = stable_inv(P,Ts,plotflag)
% Returns the stable inverse of transfer functions with only one unstable
% zero
% P: DT transfer func
%
% last update: 2013-09-07
% ============================================================
%   Copyright (c) 2008-, Xu Chen, xuchen@cal.berkeley.edu
%   Author(s): Xu Chen
% ============================================================
if nargin < 6
    P_full = P;
end
if nargin < 5
    unstable_gard = 1;
end
if nargin < 4
    z_replace = 0.6;
end
if nargin < 3
    plotflag = 'noplot';
end

[numC,denC] = tfdata(P,'v');
[zo,po,k] = tf2zp(numC,denC);
Pinv.zero = zo;
Pinv.pole = po;
Pinv.gain = k;
for ii = 1:length(zo)
    if abs(zo(ii)) >= unstable_gard
        zo_unstable = zo(ii);
    else
        if exist('den_Pinv','var')
            den_Pinv = (conv([1 -zo(ii)],den_Pinv));% (z-zo_unstable)(z-zo_stable)
        else
            den_Pinv = [1 -zo(ii)];
        end
    end
end
if exist('zo_unstable','var') && ~isempty(zo_unstable)
    Pinv.unstableZero = zo_unstable;
    num_Pinv = denC;
    for ii = 1:length(numC)
        if numC(ii) ~= 0
            gain_factor = numC(ii);
            break
        end
    end
    % the stable inverse replaces the unstable zero with a stable one
    Pinv.tf = (1+z_replace)/(1+abs(zo_unstable))/k*...
        stab_pz_cancel(...
        tf(...
        num_Pinv,...
        conv(  den_Pinv, [1 z_replace 0]  ),...
        Ts)...
        );
else
    z = tf('z',Ts);
    Pinv.tf = stab_pz_cancel(z^-(length(po)-length(zo))/P_full);
%     Pinv.tf = 1/k*...
%         stab_pz_cancel(...
%         tf(...
%         denC,...
%         den_Pinv,...
%         Ts)...
%         );
end
% G_ZPET.delaySteps = length(denC)-length(numC)+1;
try % 2011-09-18
    chkPt_Hz        = 3000;
    if 0
        chkPhase        = -angle(...
            freqresp(P_full,chkPt_Hz*2*pi)...
            ) + ...
            angle(...
            freqresp(1/Pinv.tf,chkPt_Hz*2*pi)...
            );
    else
        [mag_Pfull,pha_Pfull] = bode(P_full,chkPt_Hz*2*pi);
        [mag_Pn,pha_Pn] = bode(1/Pinv.tf,chkPt_Hz*2*pi);
        
        %     pha_Pfull = mod(pha_Pfull+180,360)-180;
        %     pha_Pn = mod(pha_Pn+180,360)-180;
        
        chkPhase_deg        = -pha_Pfull + ...
            pha_Pn;
        chkPhase_deg = mod(chkPhase_deg,360);
    end
    if 0
        %%
        figure, bode(P_full,1/Pinv.tf)
        legend('1','2')
        bode(P_full,chkPt_Hz*2*pi)
        bode(1/Pinv.tf,chkPt_Hz*2*pi)
    end
    sgn         = sign(chkPhase_deg);
    unit_deg    = Ts*chkPt_Hz*2*180; %Ts*freq_Hz*2*pi*180/pi
    delay_step  = sgn*round(chkPhase_deg/unit_deg);
    Pinv.delaySteps = max(0, delay_step);
catch
    Pinv.delaySteps = 1;
end
% C_feedforward = C_ZPET*C_Efilter;

if strncmpi(plotflag,'plot',4)
    z = tf('z',Ts);
    bode_opt = bodeoptions;
    bode_opt.FreqUnits = 'Hz';
    bode_opt.Xlim = [10 0.5/Ts];
    bode_opt.PhaseWrapping = 'On';
    bode_opt.Grid = 'On';
    
    figure;pzmap(P);title('plant');
    figure;pzmap(Pinv.tf);title('stable inverse');
    figure,bode(z^-Pinv.delaySteps/Pinv.tf,P,P*Pinv.tf,bode_opt)
    if 0
    figure,bode(z^-Pinv.delaySteps/Pinv.tf,P,P/(z^-Pinv.delaySteps/Pinv.tf),bode_opt)
    end
    legend('1/inv','P','P*inv')
end