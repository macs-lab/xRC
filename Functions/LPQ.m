function varargout = LPQ(BW_Hz,Ts,nOrder)
% function sysd_Q = LPQ(tao,Ts,nOrder)
% low pass filter design.
% ============================================================
%   Copyright (c) 2008-, Xu Chen, xuchen@cal.berkeley.edu
%   Author(s): Xu Chen
% ============================================================

% Q filter generation
% tao: time constant
% nOrder: order of Q filter
% Ts: sampling time for c2d

% tao = 0.0005;%2;%0.02;%0.0008;
% Ts  = 1/26400;
tao = 1/BW_Hz;
switch nOrder
    case 1,
        numq = 1;
        denq = [tao 1];
    case 2,
        numq = [2*tao 1];
        %         numq = [1];
        denq = [tao^2 2*tao 1];
    case 3,
        if 1
            numq = [3*tao 1];
        else
            numq = [3*tao^2 3*tao 1];
        end
        %         numq = [1];
        % numq = [tao 1];
        denq = [tao^3 3*tao^2 3*tao 1];
    case 4,
        numq = [7.7*tao^2 5.1*tao 1];
        denq = [tao^4 4.7*tao^3 7.7*tao^2 5.1*tao 1];
    case 5,
        numq = [4.7*tao^3 7.7*tao^2 5.1*tao 1];
        denq = [tao^4 4.7*tao^3 7.7*tao^2 5.1*tao 1];
    case 6,
        numq = [10*tao^2 5*tao 1];
        denq = [tao^5 5*tao^4 10*tao^3 10*tao^2 5*tao 1];
    case 7,
        numq = [10*tao^3 10*tao^2 5*tao 1];
        denq = [tao^5 5*tao^4 10*tao^3 10*tao^2 5*tao 1];
    case 8,
        numq = [5*tao^4 10*tao^3 10*tao^2 5*tao 1];
        denq = [tao^5 5*tao^4 10*tao^3 10*tao^2 5*tao 1];
end

sysc_Q = tf(numq,denq);
sysd_Q = c2d(sysc_Q,Ts,'tustin');
[numd,dend] = tfdata(sysd_Q,'v');
% numd=conv(conv([1 1],[1 1]),[6*tao/Ts+1, 1-6*tao/Ts]);
% dend=conv(conv([2*tao/Ts+1,1-2*tao/Ts],[2*tao/Ts+1,1-2*tao/Ts]),[2*tao/Ts+1,1-2*tao/Ts]);
% sysdq=tf(numd,dend,1/26400);
sysdq1 = tf(dend-numd,dend,1/26400);

if nargout == 0
    figure;
    bodeplot(tf(numq,denq));
    grid on;
    hold on;
    h=bodeplot(sysd_Q,'r--');
    setoptions(h,'FreqUnits','Hz');
    h=bodeplot(sysdq1,'k');
    legend('Q(s)','Q(z)','1-Q(z)')
    hold off;
    
%     figure;
%     pzmap(sysd_Q)
    
    %(1-Q(z))^-1
    num1q = dend;
    den1q = dend-numd;
    
    %Q*Gn^-1
    k0 = 3.744881889763779e+009;
    numqg = 4/k0/Ts/Ts*conv([1 -2 1],[6*tao/Ts+1, 1-6*tao/Ts]);
    denqg = dend;
else
    varargout{1} = sysd_Q;
end
