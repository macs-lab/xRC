function y = singen(w1,amp,phi,N,Ts)
% generates sinusoidal sequence
% w1: frequency, const
% amp: amplitude of the signal, should be a vector
% y(t) = amp(t)*sin(w1*t+phi(t))
% ============================================================
%   Copyright (c) 2008-, Xu Chen, xuchen@cal.berkeley.edu
%   Author(s): Xu Chen 
% ============================================================



for ii = (1:N)
    y(ii) = amp(ii)*sin(w1*ii*Ts+phi(ii));
end
