function BW_Hz = alpha2bw(alpha,Ts)
% function bw_Hz = alpha2bw(alpha,Ts)
% Xu Chen
% xuchen@cal.berkeley.edu
% 2011-09-22


% BW_nom      = 2*pi*BW_Hz*Ts;
% 
% alpha       = sqrt((1-tan(BW_nom/2))./(1+tan(BW_nom/2)));

BW_nom = 2*(1-alpha.^2)./(alpha^2+1);
BW_Hz = BW_nom/2/pi/Ts;