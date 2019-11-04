function varargout = zeroPhaseLowPass(Ts,varargin)
% Generates a zero phase low pass filter
%   function [Q,FFsteps,Q_delayed] = zeroPhaseQfilter(Ts,varargin)
%
%   number of preview steps required: length(zero(Q))-length(pole(Q))
%
%   Q = (a+z^-1)*(conj(a)+z^-1)/(a+1)/(conj(a)+1)*Q_LP;
%   Q = (a+z^-1)*(conj(a)+z^-1)/(a+1)/(conj(a)+1)*Q_LP;
%   Q = Q_LP^a;
%
%   It may be beneficial to locate some notches at certain
%   (e.g., resonance) frequencies.
%   In that case, choose b = cos(2*pi*Ts*freq_Hz)
% ============================================================
%   Copyright (c) 2008-, Xu Chen, xuchen@cal.berkeley.edu
%   Author(s): Xu Chen
% ============================================================

z       = tf('z',Ts);
Q       = (z^-1+2+z)/4;     % basic low pass

nFreq  = size(varargin,2)/2;

% figure;bode((1-2*c*z^-1+z^-2)*(1-2*c*z+z^2)/(2-2*c)^2)
% figure;pzmap((1-2*c*z^-1+z^-2)*(1-2*c*z+z^2)/(2-2*c)^2)
for ii = 1:nFreq
    b = varargin{(ii-1)*2 + 1};
    alpha = varargin{(ii-1)*2 + 2};
    
    % b = cos(omega) = cos(2*pi*freq_Hz*Ts)
    Q = (1  -2*alpha*b*z^-1   +alpha^2*z^-2)*...
        (1  -2*alpha*b*z      +alpha^2*z^2)/...
        (alpha^2    +1  -2*alpha*b)^2*...
        Q;
end

FFsteps     = length(zero(Q))-length(pole(Q));
Q_delayed   = z^(-FFsteps)*Q;
bw          = bandwidth(Q)/2/pi;
disp(['Bandwidth: ',num2str(bw),'[Hz]'])

if nargout == 0
    figure;     pzmap(Q);
    figure;     h = bodeplot(Q);
    setoptions(h,'FreqUnits','Hz');
    setoptions(h,'PhaseWrap','On');
    xlim([10,0.5/Ts])
    title('zero phase low pass Q filter')
elseif nargout == 3
    varargout{1} = Q;
    varargout{2} = FFsteps;
    varargout{3} = Q_delayed;
else
    error('Wrong number of outputs')
end