function y = signalgen_ti(Ts,w,A,Pha,n,N,m,plotflag,noiseflag)
% function y = signalgen_ti(Ts,w,A,Pha,n,N,m,plotflag,noiseflag)
% 
% Ts; sampling time
% w: sin(w*Ts*k)
% n: number of narrow band disturbances
% N: signal length
% ti: time invariant mag and phase
% A: Amplitude vector
% Pha: Phase vector
% noiseflag: flag for adding noise to the amplitude and phases

% ============================================================
%   Copyright (c) 2008-, Xu Chen, xuchen@cal.berkeley.edu
%   Author(s): Xu Chen 
% ============================================================



% Ts = 1/26400;
% Ts = 60/26400;
% N = 1/60*Ts;
t       = N*Ts;

% w1 = 0.8*pi/Ts;
% w1 = 0.08*pi/Ts;
% w2 = 0.05*pi/Ts;
% w2 = 0.2*pi/Ts;
% f1 = w1/2/pi;%frequency in Hz
% f2 = w2/2/pi;
%
N       = N+1;
wN      = zeros(n,1);
wperiod = zeros(n,1);
amp     = zeros(N,n);
phi     = zeros(N,n);
x       = zeros(N,n);
yperiod = 1;
for ii = 1:n
    [wN(ii), wperiod(ii)]   = periodCal(w(ii)*Ts, Ts, 0.001);
    if strcmp(noiseflag,'noiseOn')
        amp(220*m+1:end,ii) = A(ii)*ones(N-220*m,1)+rand(N-220*m,1)*1e-3;
        %                               first m revolutions no disturbance
        phi(220*m+1:end,ii) = Pha(ii)*ones(N-220*m,1)+rand(N-220*m,1);
    else
        amp(220*m+1:end,ii) = A(ii)*ones(N-220*m,1);
        phi(220*m+1:end,ii) = Pha(ii)*ones(N-220*m,1);
    end
    % amp(k)*sin(w1*Ts*k+phi(k)) k=1:N
    x(:,ii) = singen(w(ii),amp(:,ii),phi(:,ii),N,Ts);
    yperiod = lcm(yperiod,wN(ii));
end

y = sum(x,2);

if strncmpi(plotflag,'plot',4)
    for ii = 1:n
        figure;
        plot([0:N-1]'*Ts,x(:,ii));title(sprintf('The %d-th narrow band component',ii));
        xlabel(sprintf('time [sec], T = %f', wperiod(ii)));
        ylabel 'Amplitude';
    end
    figure;
    plot([0:N-1]'*Ts,y);title 'The multiple narrow band disturbance';
    xlabel(sprintf('time [sec], T = %f', yperiod*Ts));
    ylabel 'x';
end
