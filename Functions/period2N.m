function delay_step = period2N(Ts,freq_Hz)
% function delay_step = period2N(Ts,freq_Hz)
%   freq_Hz: freq in Hz
% Author: Xu Chen xuchen@cal.berkeley.edu
% Initial Version: 2011-06-24
unit_deg    = Ts*freq_Hz*2*180; %Ts*freq_Hz*2*pi*180/pi
pha = 2*pi;
req_deg     = pha*180/pi;
delay_step  = round(req_deg./unit_deg);
