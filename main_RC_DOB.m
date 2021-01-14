% =================================================
% RC_DOB: repetitive control via the DOB structure
% algorithms:   youla equivalence
%               IMP
%               diophantine eq
%
% Ref:
%       RC_DOB.m
% =================================================
% 2011-09-18 Xu Chen chx@uw.edu
% ============================================================
%   Copyright (c) 2008-, Xu Chen, chx@uw.edu
%   Author(s): Xu Chen
% ============================================================
folder = fileparts(which(mfilename));
addpath(genpath(folder));

dataLoad;
SW_SIMPLE = 1;
revo_transient = 5*60/7200;

%% User input
if SW_SIMPLE
    g_nb = 0;
else
    flag_NB = input('NB dist on? [yes-enter; no-any input]:');
    % switch_DOB = isempty(flag_DOB)? 1:0;
    if ~isempty(flag_NB)
        g_nb = 0;
    else
        g_nb = 1;
    end
end

if SW_SIMPLE
    flag_NB_method = 1;
else
    flag_NB_method = input(...
        'NB generation method? [sin-enter; white noise excitation-any input]:');
end

flag_Dist = input(...
    'Other disturbance on? [yes-enter; no-any input]:');
if ~isempty(flag_Dist)
    g_torque  = 0;
    g_sensor  = 0;
    g_flutter = 0;
    g_rro     = 1;
else
    g_torque  = 1;
    g_sensor  = 1;
    g_flutter = 1; % disk flutter should be excluded here since it's the
    % main source of narrow band disturbances and varies
    % between tracks
    %     g_rro     = 1;
    g_rro     = 1; % in reality, rro is always compensated, e.g., by AFC
end

flag_DOB = input('DOB on? [yes-enter; no-any input]:');
if ~isempty(flag_DOB)
    switch_DOB = 0;
else
    switch_DOB = 1;
end

flag_Q = input(...
    'Peak filter form: [1-enter; 2(perfect Q for 1 delay)-any input]:');

flag_SW_SAVE_DATA = input(...
    'SAVE DATA?: [no-enter; yes-any input]:');
if isempty(flag_SW_SAVE_DATA)
    SW_SAVE_DATA = 0;
else
    SW_SAVE_DATA = 1;
end

adap_method = 0;
SimModel = 'RC_DOB_mdl';

flag_DOBperform = input(...
    'Analyze DOB performance? [No-enter; Yes-any input]:');

flag_Qcode = input(...
    'Adapt shaping coefficient in Q? [No-any input; Yes-enter]:');
if isempty(flag_Qcode)
    filterQ_gain = 0;
    codeQ_gain = 1;
else
    filterQ_gain = 1;
    codeQ_gain = 0;
end
distDefine;
%% Feedback controllers
FBcontroller;
%% Simulation
% Option for figure plot
close_figure = 0; % (1: close figure, 0: keep all figures)

% Loop gain perturbation
loopgain_pt = PlantData.DeltaKp;

% For robust switch
if (robust_sw==1)
    gnum         = length(loopgain_pt);
    [nn,mm,pnum] = size(PlantData.Pfpert);
elseif (robust_sw==2)
    gnum = 1;
    pnum = 1;
else
    gnum = length(loopgain_pt);
    pnum = 1;
end
% for debug
% if pnum > 3
%    pnum = 3;
% end
num_servo = PlantData.num_servo;

% Sampling time for plant input
% Tu = PlantData.Tu;
Tu = Ts;

% Simulation Time
if ~exist('Tsim','var')
    Tsim = NmaxTs*Ts;
end
%%
% G_nominal
Pn = PlantData.Pn;
Ts = PlantData.Ts;
[num_Pn, den_Pn] = ss2tf(Pn.a, Pn.b, Pn.c, Pn.d);
Pn_tf = tf(num_Pn, den_Pn);
% figure;bode(Pn,Pn_tf);legend('Pn','Pn tf');title('tf loses the delay')
% 1.137e-013 s + 3.745e009
% -------------------------
% s^2 + 565.5 s + 3.198e005
Pn_d = c2d(Pn,Ts,'zoh');
[num_Pn_d, den_Pn_d] = ss2tf(Pn_d.a, Pn_d.b, Pn_d.c, Pn_d.d);
Pn_d_tf = tf(num_Pn_d, den_Pn_d,Ts);
% figure;bode(Pn,Pn_d_tf);legend('Pn','Pn d tf')
% figure;bode(Pn_d_tf,Pn_tf);legend('Pn d tf','Pn tf')
% 1.448 z^2 + 3.685 z + 0.1836
% ----------------------------
%  z^3 - 1.978 z^2 + 0.9788 z
%
% Sampling time: 3.7879e-005

if 1 % RC DOB configuration 2011-11-20
    if isempty(flag_Qcode)
        
        alpha = input(' alpha = (0.999 by default)');%0.999;
        if isempty(alpha)
            alpha = 0.999;
        end
        % simulink initialization
        alpha_init = 1;
        alpha_rate = 0.9;
        alpha_rate1 = 0.9;
        alpha_rate2 = 0.9;
%         alpha_end = 0.999; % 2012-01-25
        alpha_end = alpha;
        alpha_judge = 0.01;
        alpha_middle = alpha_end;
    else
        alpha_init = input('    alpha_init = (1 by default)');%1;%0.999;
        if isempty(alpha_init)
            alpha_init = 1;
        end
        
        alpha_middle = input('    alpha_middle = (0.99 by default)');%0.999;
        if isempty(alpha_middle)
            alpha_middle = 0.99;
        end
        
        alpha_end = input('    alpha_end = (0.999 by default)');%0.999;
        if isempty(alpha_end)
            alpha_end = 0.999;
        end
        %           0.99 gives very good transient, but amplification of around
        %               1800, needs a reduced bw ZPET low-pass filter
        %           <0.99 gives not good result, has to reduce bw of ZPET lp
        %
        
        alpha_judge = input('    alpha_converge_judge_ratio = (0.01 by default)');%0.999;
        if isempty(alpha_judge)
            alpha_judge = 0.01;
        end
        
        alpha_rate1 = input('    alpha_rate_stage1 = (0.9 by default)');%0.9;
        if isempty(alpha_rate1)
            alpha_rate1 = 0.9;
        end
        
        alpha_rate2 = input('    alpha_rate_stage2 = (0.95 by default)');%0.9;
        if isempty(alpha_rate2)
            alpha_rate2 = 0.95;
        end
        alpha = 0.999;
    end
else
    alpha_init = 0;
    alpha_end = 0;
    alpha_rate = 0.9;
    alpha = 0;
end
%%
lowP = LPQ(1000*2*pi,Ts,4);
if 0
    %%
    for ii = 3:8
        loww(:,:,ii-2) = LPQ(1000*2*pi,Ts,ii)
        lowww{ii-2} = LPQ(1000*2*pi,Ts,ii)
    end
    figure, bodeplot(loww,bodeopt)
    figure, xbodeplot(lowww)
end
% Run simulation
% close all
for kkk = 1:2
    if ~exist('nosim')
        clear NRPE3sigma RPEpp Result;
        
        idx = 1;
        for jj=1:gnum
            for ii=1:pnum
                loopgain = 1+loopgain_pt(jj);
                %             fprintf('Model No.%d, loopgain = %4.2f : ',ii,loopgain);
                sys_P = PlantData.Pfpert(1,1,ii);
                sysd_P = c2d(ssbal(sys_P),Tu); % discretize with input delay (Tc -> Tu, Ver.3.1)
                %% Nominal Plant
                %             PnZPETinverse = NominalPselect(Pn,sysd_P,Ts,'zpet','plot');
                PnZPETinverse = NominalPselect(Pn,sysd_P*sys_C2,Ts,'sInv','plot');
                %             PnZPETinverse = NominalPselect(Pn,sysd_P,Ts,'sInv','plot');
                G_n_inv0 = PnZPETinverse.tf;
                d_DOB = PnZPETinverse.delaySteps;
                m = d_DOB;
                %% Q
                %                 [Q,Sadd] = RC_DOB(120,Ts,1);
%                 baseFreq_Hz,Ts,m, alpha, SW_LPLP, SW_TEST, pair_LPLP
                if alpha < 0.9
%                     function varargout = RC_DOB(baseFreq_Hz,Ts,m, alpha,
%                     SW_LPLP, SW_TEST, pair_LPLP)
                    [Q,Sadd,LPLP,n_lplp] = RC_DOB(120,Ts,m,alpha_end,1,0,4);
                else
                    [Q,Sadd,LPLP,n_lplp] = RC_DOB(120,Ts,m,alpha_end);
                end
                %%
                tic;
                if kkk == 1
                    switch_DOB = 0;
                else
                    switch_DOB = 1;
                end
                sim(SimModel,Tsim);
                etime = toc;
                idxt = find(PESDataT.time >= 10e-3);        
                PES  = PESDataT.signals.values(idxt)*1e2;	% [%TP]
                if close_figure, close, end
                [rpe,nrpe,val] = PES_PlotTD(PES,num_servo,'%TP');
                drawnow
                %             figsize(500,320,'keep'); % comment out (Ver.3.1)
                val_nrpe3sigma = val.NRPE6sigma/2;
                val_rpepp      = val.RPEpp;
                fprintf('NRPE 3sigma = %f, RPEpp = %f',val_nrpe3sigma,val_rpepp);
                fprintf(' (elapsed time = %g (s))', etime);
                fprintf('\n')
                Result(idx,:)   = [ii,loopgain,val_nrpe3sigma,val_rpepp];
                NRPE3sigma(idx) = val_nrpe3sigma;
                RPEpp(idx)      = val_rpepp;
                idx = idx + 1;
                if kkk == 1
                    TMR_nocomp = getTMR(PESDataT.signals.values((t_Qon+revo_transient)/Ts:end));
                    PES_nocomp = PESDataT;
                    specPES_nocomp=specCal(PESDataT.signals.values((t_Qon+revo_transient)/Ts:end),1/Ts);
                end
            end
        end
    end
end
figure,
xbodeplot({sysd_P*sys_C2,z^-m/PnZPETinverse.tf})
legend('Plant with notch filters', 'Nominal model','Location','Best');

if SW_SAVE_DATA
    fig_name=['bode_P_Pn'];
    hgsave(['Src/RC_DOB/',fig_name])
    saveas(gcf,['Src/RC_DOB/',fig_name,'.emf'])
end

figure,
xbodeplot({Q})
title 'Q filter'
if SW_SAVE_DATA
    fig_name=['bode_Q'];
    hgsave(['Src/RC_DOB/',fig_name])
    saveas(gcf,['Src/RC_DOB/',fig_name,'.emf'])
end

if gnum*pnum > 1
    % Plot summary
    NRPE3sigma2 = reshape(NRPE3sigma,pnum,gnum);
    RPEpp2      = reshape(RPEpp,pnum,gnum);
    idxp        = 1:pnum;
    
    % NRPE3sigma
    figure;
    plot(idxp,NRPE3sigma2(:,1),'r*-',...
        idxp,NRPE3sigma2(:,2),'go-',...
        idxp,NRPE3sigma2(:,3),'bs-')
    set(gca, 'XTick', [1:pnum])
    set(gca, 'YLim', [7 10])
    xlabel('Model number'), ylabel('NRPE 3sigma [percent/track]')
    legend('Nominal loopgain','loopgain = +10%','loopgain = -10%')
    title('None Repeatable Position Error (NRPE)')
    
    % RPEpp
    figure;
    plot(idxp,RPEpp2(:,1),'r*-',...
        idxp,RPEpp2(:,2),'go-',...
        idxp,RPEpp2(:,3),'bs-')
    set(gca, 'XTick', [1:pnum])
    set(gca, 'YLim', [10 14])
    xlabel('Model number'), ylabel('RPEpp [percent/track]')
    legend('Nominal loopgain','loopgain = +10%','loopgain = -10%')
    title('Repeatable Position Error Peak-to-Peak (RPEpp)')
end
%% Disturbance rejection result
ylim_max = max(max(abs([PESDataT.signals.values;PES_nocomp.signals.values])))*100;
ylim_bd = [-ylim_max,ylim_max];

TMR = getTMR(PESDataT.signals.values((t_Qon+revo_transient)/Ts:end));

figure;plot(PESDataT.time*120,100*PESDataT.signals.values);grid;
xlabel('Revolution');ylabel('PES (%TP)');
ylim(ylim_bd)
text(t_Qon*120,min(100*PESDataT.signals.values)+2,'\leftarrow Compensation starts',...
    'FontSize',12,'color','r');
line([t_Qon*120 t_Qon*120],ylim_bd,'Color','r','LineStyle','-','LineWidth',3)
if SW_SAVE_DATA
    fig_name=['t_PES_w_comp'];
    hgsave(['Src/RC_DOB/',fig_name])
    saveas(gcf,['Src/RC_DOB/',fig_name,'.emf'])
end

figure;
subplot(211)
plot(...
    PES_nocomp.time*120,100*PES_nocomp.signals.values);grid;
ylabel('PES (%TP)');
ylim(ylim_bd)
legend('w/o compensation','location','best')
subplot(212)
plot(...
    PESDataT.time*120,100*PESDataT.signals.values);grid;
ylim(ylim_bd)
xlabel('Revolution');ylabel('PES (%TP)');
text(t_Qon*120,min(100*PESDataT.signals.values)+2,'\leftarrow Compensation starts',...
    'FontSize',12,'color','r');
line([t_Qon*120 t_Qon*120],ylim_bd,'Color','r','LineStyle','-','LineWidth',3)
legend('w/ compensation','location','best')

if SW_SAVE_DATA
    fig_name=['t_PES_subplot_copmpare'];
    hgsave(['Src/RC_DOB/',fig_name])
    saveas(gcf,['Src/RC_DOB/',fig_name,'.emf'])
end

specPES=specCal(PESDataT.signals.values((t_Qon+revo_transient)/Ts:end),1/Ts);
figure;
plot(specPES.f,specPES.amp)
% title('PES Amplitude Spectrum after compensation')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
xlim([1,0.5/Ts]);
legend(...
    sprintf('w/ compensation 3\\sigma = %.2f %%TP',TMR*100))
if SW_SAVE_DATA
    fig_name=['spec_PES_w_comp'];
    hgsave(['Src/RC_DOB/',fig_name])
    saveas(gcf,['Src/RC_DOB/',fig_name,'.emf'])
end
% ylim([0,0.06]);
%
figure;
plot(specPES_nocomp.f,specPES_nocomp.amp)
% title('PES Amplitude Spectrum after compensation')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
xlim([1,0.5/Ts]);
legend(...
    sprintf('w/o compensation 3\\sigma = %.2f %%TP',TMR_nocomp*100))

if SW_SAVE_DATA
    fig_name=['spec_PES_wo_comp'];
    hgsave(['Src/RC_DOB/',fig_name])
    saveas(gcf,['Src/RC_DOB/',fig_name,'.emf'])
end

figure;
plot(specPES.f,specPES.amp)
hold on
plot(specPES_nocomp.f,specPES_nocomp.amp,'r:')
title('PES Amplitude Spectrum')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
legend(sprintf('w/ compensation 3\\sigma = %.2f %%TP',TMR*100),...
    sprintf('w/o compensation 3\\sigma = %.2f %%TP',TMR_nocomp*100))
xlim([1,0.5/Ts]);

if SW_SAVE_DATA
    fig_name=['spec_PES_copmpare'];
    hgsave(['Src/RC_DOB/',fig_name])
    saveas(gcf,['Src/RC_DOB/',fig_name,'.emf'])
end

figure,
subplot(211)
plot(specPES_nocomp.f,specPES_nocomp.amp)
yl = ylim;
xlim([1,0.5/Ts])
legend(...
    sprintf('w/o compensation 3\\sigma = %.2f %%TP',TMR_nocomp*100))
ylabel 'Magnitude'
subplot(212)
plot(specPES.f,specPES.amp)
legend('w/ compensation')
ylabel 'Magnitude'
xlabel 'Frequency (Hz)'
ylim(yl)
xlim([1,0.5/Ts])
legend(sprintf('w/ compensation 3\\sigma = %.2f %%TP',TMR*100)...
    )

if SW_SAVE_DATA
    fig_name=['spec_PES_copmpare_lin_subplot_full'];
    hgsave(['Src/RC_DOB/',fig_name])
    saveas(gcf,['Src/RC_DOB/',fig_name,'.emf'])
end

figure,
subplot(211)
plot(specPES_nocomp.f,specPES_nocomp.amp)
yl = ylim;
xlim([1,2000])
% legend('w/o compensation')
legend(...
    sprintf('w/o compensation 3\\sigma = %.2f %%TP',TMR_nocomp*100))
ylabel 'Magnitude'
subplot(212)
plot(specPES.f,specPES.amp)
ylabel 'Magnitude'
% legend('w/ compensation')
legend(sprintf('w/ compensation 3\\sigma = %.2f %%TP',TMR*100)...
    )
xlabel 'Frequency (Hz)'
ylim(yl)
xlim([1,2000])

if SW_SAVE_DATA
    fig_name=['spec_PES_copmpare_lin_subplot_main'];
    hgsave(['Src/RC_DOB/',fig_name])
    saveas(gcf,['Src/RC_DOB/',fig_name,'.emf'])
end

figure;
semilogx(specPES.f,specPES.amp)
hold on
semilogx(specPES_nocomp.f,specPES_nocomp.amp,'r:')
title('PES Amplitude Spectrum')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
legend(sprintf('w/ compensation 3\\sigma = %.2f %%TP',TMR*100),...
    sprintf('w/o compensation 3\\sigma = %.2f %%TP',TMR_nocomp*100))
xlim([10,0.5/Ts]);


if SW_SAVE_DATA
    fig_name=['spec_PES_copmpare_log_full'];
    hgsave(['Src/RC_DOB/',fig_name])
    saveas(gcf,['Src/RC_DOB/',fig_name,'.emf'])
end

if 0
    specDhat=specCal(d_hat.signals.values,1/Ts);
    figure;
    semilogx(specDhat.f,specDhat.amp)
    title('DOB Dist Estimate Amplitude Spectrum')
    xlabel('Frequency [Hz]')
    ylabel('Signal Power')
    xlim([1,2e4]);
end
%% DOB performance
if ~isempty(flag_DOBperform)
    if 1
        %%
        DOBperform(sys_C1,1,sysd_P*sys_C2,Pn_d_tf,Q,G_n_inv0,d_DOB,z);
        if 0 % robust stability test
            %%
            alpha_end = 0.999;
            SW_LPLP = 1;
            [Q,Sadd,LPLP,n_lplp] = RC_DOB(120,Ts,m,alpha_end,SW_LPLP);
            [S_wDOB_lplp,S_woDOB_lplp] = DOBperform(...
                sys_C1,1,sysd_P*sys_C2,Pn_d_tf,Q,G_n_inv0,d_DOB,z);
            SW_LPLP = 0;
            [Q,Sadd,LPLP,n_lplp] = RC_DOB(120,Ts,m,alpha_end,SW_LPLP);
            [S_wDOB,S_woDOB] = DOBperform(...
                sys_C1,1,sysd_P*sys_C2,Pn_d_tf,Q,G_n_inv0,d_DOB,z);
            %% try a single plot first
            T_wDOB_lplp = 1-S_wDOB_lplp;
            T_woDOB_lplp = 1-S_woDOB_lplp;
            T_wDOB = 1-S_wDOB; 
            T_woDOB = 1-S_wDOB;
            
            figure, subplot(211)
            xbodemag({S_wDOB_lplp,S_woDOB_lplp})
            xlim([10,0.5/Ts]); legend('w/ RDOB','w/o RDOB','loction','best')
            title('Sensitivity functions')
            subplot(212)
            xbodemag({1/T_wDOB_lplp,1/T_woDOB_lplp,1/T_woDOB},1.5,'linestyleOn',0.5/Ts,2000)
            xlim([10,0.5/Ts]);
            title('Uncertainty bound for robust stability')
            legend('w/ RDOB','w/o RDOB','w/ RDOB; no lowpass in Q','loction','best')
            %% separate plot 1
            figure
            xbodemag({S_wDOB_lplp,S_woDOB_lplp})
            xlim([10,0.5/Ts]); legend('w/ RDOB','w/o RDOB','loction','best')
%             title('Sensitivity functions')
            %% separate plot 2
            figure,
            subplot(211)
            xbodemag({1/T_wDOB_lplp,1/T_woDOB_lplp},1.5,'linestyleOn',0.5/Ts,2000)
            xlim([10,0.5/Ts]);
            ylim([-32,50]);
            grid off;
            xlabel '';
%             title('Uncertainty bound for robust stability')
            legend('w/ RDOB','w/o RDOB','loction','best')
            subplot(212)
            xbodemag({1/T_woDOB,1/T_woDOB_lplp},1.5,'linestyleOn',0.5/Ts,2000)
            xlim([10,0.5/Ts]);
            ylim([-32,50]);
            grid off;
%             title('Uncertainty bound for robust stability')
            legend('w/ RDOB; but no lowpass filter in Q','w/o RDOB','loction','best')
        end
    end
end
%% SPR test for the IMP
figure,
plot(Q_shaping_coef.signals.values)
xlabel 'Sample'
ylabel '\alpha'

if SW_SAVE_DATA
    fig_name=['t_Q_shaping_coef'];
    hgsave(['Src/RC_DOB/',fig_name])
    saveas(gcf,['Src/RC_DOB/',fig_name,'.emf'])
end

figure,
plot(Q_shaping_coef.time*120,Q_shaping_coef.signals.values)
xlabel 'Revolution'
ylabel '\alpha'

if SW_SAVE_DATA
    fig_name=['t_Q_shaping_coef_revo'];
    hgsave(['Src/RC_DOB/',fig_name])
    saveas(gcf,['Src/RC_DOB/',fig_name,'.emf'])
end
