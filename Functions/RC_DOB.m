function varargout = RC_DOB(baseFreq_Hz,Ts,m, alpha, SW_LPLP, SW_TEST, pair_LPLP)
% function varargout = RC_DOB(baseFreq_Hz,Ts,m, alpha, SW_LPLP, SW_TEST,
% pair_LPLP)
% 
% 2011-08-30
% Xu Chen
% xuchen@cal.berkeley.edu
% repetitive control using a DOB
if nargin < 7
    pair_LPLP = 2;
end
if nargin < 6
    SW_TEST = 0;
end
if nargin < 5
    SW_LPLP = 1;
end
if nargin < 4
    alpha = 0.999;
end
if nargin < 3
    m = 1;
end
if nargin < 2
    Ts = 1/26400;
end
if nargin < 1
    %     N = 220;
    baseFreq_Hz = 120;
end
% NBw = 1500*2*pi;
% SW_LPLP = 1;
NBw = baseFreq_Hz*2*pi;

% m = 1;
% Ts = 1/26400;
% if nargin < 1
%     N = 220;
N = period2N(Ts,NBw/2/pi)-1; % for some reason I used minus one
% end

z = tf('z',Ts);
if ~SW_TEST
    %% narrow-band case
    % alpha = 0.999;
    % lattice_bw = 50;
    lattice_bw = alpha2bw(alpha,Ts);
    lbtrue = 2*cos(NBw*Ts);
    % Aq = 1-2*cos(2*pi*w*Ts)*alpha*z^-1+ alpha^2*z^-2;
    % Bq =
    
    z       = tf            ('z',Ts);
    if 1
        [G,Q]  = lattice_filter(NBw/2/pi, lattice_bw, Ts);
        %     Q      = Xu_optQ_dirCal(Q_fc,Q1,m,Ts); % aop.optQ_NdelayCompensate
        %     Q = Q/2;
        %     Q = tf(1,1,Ts);
        %     Q = LPQ(2000,Ts,3);
    else
        NF   = NotchProd     (alpha,1,lbtrue/2,Ts,'noplot');
        Q    = minreal       ((1-NF)*z);
    end
    Sadd = 1-z^-m*z^-N*Q;
    % Sadd = 1-z^-m*z^-N*Q;
    bodeopt = fun_defn_bode_opt(10,Ts);
    if nargout == 0
        figure, xbodemag({Sadd,Q})
        grid
        legend('1-z^{-m}Q(z^{-1})','Q(z^{-1})','location','best')
        
        figure, xbodemag({Sadd})
        grid, title '1-z^{-m}Q(z^{-1})'
        
        figure, pzplot(Sadd), title '1-z^{-m}Q(z^{-1})'
    end
    % reasons: 1/(1-z^-N) is the internal model for a group of periodic
    % disturbances (see repetitive control). Q is 1 at its center frequency,
    % and decaying gains on the two sides. At center frequency, 1 gives perfect
    % disturbance rejection, on the two sides, the perfect cancellation zero
    % has decreasing magnitudes.
end
%% perfect RC_DOB
% ref: optMultiQ.m
% alpha = 0.9;
% note N is N-1 in repetitive control
if SW_TEST
    N = 10;
%     SW_LPLP = 0;
end
A = [1 zeros(1,N) -1];
if 0
    Bq_fxd = [1 1];
    % Bq_fxd = [1 0 -1];
else
    Bq_fxd = 1;
end
Coef_delay = [zeros(1,m), 1];

if 1
    Aq = [1 zeros(1,N) -alpha^(N+1)];
else
    Aq = [1 -2*alpha*cos( NBw(1)*Ts ) alpha^2];
end

if 0
    [Bq_add,K,nBq_add,nK] = bezoutd(A,Coef_delay,1,Bq_fxd,Aq,1e-6);
    % sometimes gives some strange solutions
    Q = tf(conv(Bq_fxd,Bq_add),Aq,Ts,'Variable','z^-1');
else
    if SW_LPLP
        %% linear phase low-pass filter
        if 0
            Fs = 26400;  % Sampling Frequency
            
            Fpass = 2000;            % Passband Frequency
            Fstop = 3000;            % Stopband Frequency
            Dpass = 0.057501127785;  % Passband Ripple
            Dstop = 0.031622776602;  % Stopband Attenuation
            dens  = 20;              % Density Factor
            
            % Calculate the order from the parameters using FIRPMORD.
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
            
            % Calculate the coefficients using the FIRPM function.
            b  = firpm(N, Fo, Ao, W, {dens});
            Hd = dfilt.dffir(b);
            
            lplp = tf(Hd.Numerator,[1 zeros(1,length(Hd.Numerator))],1/Fs,'Variable','z^-1');
            figure, bodeplot(lplp);
        else
            %%
            %     Notch frequency:
            %     4210      4210/13200
            %     8410      8410/13200*-1
            %     12200     12200/13200*-1
            %             [LP_0phase.tf, LP_0phase.delay, LP_0phase.tf_causal] = ...
            %                 zeroPhaseLowPass(...
            %                 Ts,...
            %                 cos(12200/13200*pi),    1,...
            %                 cos(8400/13200*pi),     1,...
            %                 cos(4210/13200*pi),     1); % 0.98 does not change the bandwidth
            if pair_LPLP == 2
            [LP_0phase.tf, LP_0phase.delay, LP_0phase.tf_causal] = ...
                zeroPhaseLowPass(...
                Ts,...
                cos(12200/13200*pi),    1,...
                cos(8400/13200*pi),     1);
            elseif pair_LPLP == 3
                [LP_0phase.tf, LP_0phase.delay, LP_0phase.tf_causal] = ...
                zeroPhaseLowPass(...
                Ts,...
                cos(12200/13200*pi),    1,...
                cos(8400/13200*pi),     1,...
                cos(6400/13200*pi),     1);
            elseif pair_LPLP == 4
                [LP_0phase.tf, LP_0phase.delay, LP_0phase.tf_causal] = ...
                zeroPhaseLowPass(...
                Ts,...
                cos(12200/13200*pi),    1,...
                cos(8400/13200*pi),     1,...
                cos(6400/13200*pi),     1,...
                cos(6400/13200*pi),     1);
            else
                [LP_0phase.tf, LP_0phase.delay, LP_0phase.tf_causal] = ...
                zeroPhaseLowPass(...
                Ts,...
                cos(12200/13200*pi),    1,...
                cos(8400/13200*pi),     1);
            end
            %     %%
            %
            %     [LP_0phase.tf, LP_0phase.delay, LP_0phase.tf_causal] = ...
            %         zeroPhaseLowPass(...
            %         Ts,...
            %         cos(12200/13200*pi),    1,...
            %         cos(8400/13200*pi),     1,...
            %         cos(6210/13200*pi),     1,...
            %         cos(5210/13200*pi),     1,...
            %         cos(4210/13200*pi),     1,...
            %         cos(3010/13200*pi),     1,...
            %         cos(2010/13200*pi),     1,...
            %         cos(1010/13200*pi),     1);
            %     figure,bodeplot(1-LP_0phase.tf_causal,bode_opt)
            %
            %
            % High pass is also possible, for output DOB
            %
            %     [LP_0phase.tf, LP_0phase.delay, LP_0phase.tf_causal] = ...
            %         zeroPhaseHighPass(...
            %         Ts,0,...
            %         cos(50/13200*pi),    1);
            % figure,bodeplot(sys_C2,LP_0phase.tf_causal,bode_opt)
            lplp          = stab_pz_cancel(LP_0phase.tf_causal);
            
            %     figure, bodeplot(adap_prefilter)
            %     figure, bodeplot(adap_prefilter,bode_opt)
            %     figure, bodeplot(adap_prefilter,adap_temp,'r--',bode_opt)
        end
        if 0
            %%
            figure, bodeplot(lplp)
        end
    end
    if m == 1
        if 1
            Bq = [zeros(1,N+1-m), (1-alpha^(N+1))]; % normal
        else
            Bq = [-1, zeros(1,N-1), (1-alpha^(N+1)), 1]; % high freq constraint
        end
        %     elseif m == 2
    else
        Bq = (1-alpha^(N+1))*[zeros(1,N+1-m), 1];
    end
    if SW_LPLP
        Bq = (1-alpha^(N+1))*[zeros(1,N+1-m-order(lplp)/2), 1]; % 2011-09-19
        Q = ( tf(Bq,Aq,Ts,'Variable','z^-1')*lplp );
        %         Q = stab_pz_cancel( tf(Bq,Aq,Ts,'Variable','z^-1')*lplp );
    else
        lplp = tf(1,1,Ts);
        Q = tf(Bq,Aq,Ts,'Variable','z^-1');
    end
end
if 0
    %%
    figure, bodeplot(Q)
end
% Q = Q*LPQ(4000*2*pi,Ts,3);
IMP = tf(A,[1 zeros(1,length(A)-1)],Ts,'Variable','z^-1');
% Sadd = 1-z^-m*z^-N*Q;
Sadd = 1-z^-m*Q;

if nargout == 0
    figure, xbodemag({Sadd,Q})
    grid
    ylim([-150,5])
    legend('1-z^{-m}Q(z^{-1})','Q(z^{-1})','location','best')
    
    figure, subplot(211)
    xbodemag({Sadd})
    title '1-z^{-m}Q(z^{-1})'
    ylim([-150,5])
    subplot(212)
    xbodemag({Q})
    title 'Q(z^{-1})'
    
%     if SW_TEST
%         figure, subplot(211)
%         bodemag(Sadd)
%         title '1-z^{-m}Q(z^{-1})'
%         ylim([-150,5])
%         subplot(212)
%         bodemag(Q)
%         ylim([-150,5])
%         title 'Q(z^{-1})'
%     end
    
    figure, pzplot(Q)
    title 'Q(z^{-1})'
    
    figure, pzplot(Sadd,IMP)
    legend('1-z^{-m}Q(z^{-1})','1-z^{-N}','location','best')
    
    figure, pzplot(Sadd), title '1-z^{-m}Q(z^{-1})'
    
    figure, step(Sadd)
    figure, impulse(Sadd)
elseif nargout == 1
    varargout{1} = Q;
elseif nargout == 2
    varargout{1} = Q;
    varargout{2} = Sadd;
elseif nargout == 3
    varargout{1} = Q;
    varargout{2} = Sadd;
    varargout{3} = lplp;
elseif nargout == 4
    varargout{1} = Q;
    varargout{2} = Sadd;
    varargout{3} = lplp;
    varargout{4} = order(lplp)/2;
end
if nargout ~= 0
    figure, subplot(211)
    xbodemag({Sadd})
    title '1-z^{-m}Q(z^{-1})'
    ylim([-150,5])
    subplot(212)
    xbodemag({Q})
    title 'Q(z^{-1})'
    
    figure, 
    xbodemag({Sadd})
    title '1-z^{-m}Q(z^{-1})'
    ylim([-150,5])
    figure
    xbodemag({Q})
    title 'Q(z^{-1})'
    ylim([-30,5])
end
if 0
    %%
    close all
end