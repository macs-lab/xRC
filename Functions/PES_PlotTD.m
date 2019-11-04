function [RPE,NRPE,Val] = PES_PlotTD(PES,num_servo,unit)
%PES_PlotTD  Plot PES Time History Data
%
%   PES_PlotTD(PES,num_servo,unit) plots RPE/NRPE/TPE from PES time history
%   data
%
%   PES       : PES time history data
%   num_servo : Number of servo sectors per rotation
%   unit      : (optional) unit strings eg. '%TP'
%
%   Plot Line colors:
%   Red   : RPE, repetitive PES
%   Green : TPEmax/min
%   Blue  : NRPEmax/min
%   Purple: NRPE6sigma
%
%   Copyright (c) 2005, MSS benchmark working group
%   Ver.1.0, 2005-08-16 T. Hara

if ~exist('unit')
    unit = '';
end

rd = floor(length(PES)/num_servo);

PES2 = reshape(PES(1:num_servo*rd),num_servo, rd);
RPE = mean(PES2,2);% repetitive PES

RPE2 = reshape(RPE*ones(1,rd+1),num_servo*(rd+1),1);
RPE2 = RPE2(1:length(PES));% stack for calculating the NRPE

NRPE = PES - RPE2;% non-repetitive PES

PES4 = reshape(NRPE(1:num_servo*rd),num_servo,rd);

NRPEstd = (std(PES4,0,2));% standard deviation in the vectors
NRPEmax = (max(PES4,[],2));
NRPEmin = (min(PES4,[],2));
TPEmax = (max(PES2,[],2));
TPEmin = (min(PES2,[],2));
Val.RPEpp = max(RPE)-min(RPE);
Val.TPEpp = max(PES)-min(PES);
Val.NRPE6sigma = sqrt(mean(NRPEstd.^2))*6;

Sct=[0:num_servo]';
figure;
h_ = stairs(Sct,[RPE; RPE(1)],'r');
set(h_,'LineWidth',1.5);	% 1.5pt ~= 0.529mm
hold on
h_ = stairs(Sct,[NRPEstd; NRPEstd(1)]*6,'m');
set(h_,'Color',[0.75 0 0.75]);	% Purple
h_ = stairs(Sct,[TPEmax; TPEmax(1)],'g');
set(h_,'Color',[0 0.5 0]);	% Dark Green
h_ = stairs(Sct,[TPEmin; TPEmin(1)],'g');
set(h_,'Color',[0 0.5 0]);	% Dark Green
stairs(Sct,[NRPEmax; NRPEmax(1)],'b');
stairs(Sct,[NRPEmin; NRPEmin(1)],'b');
legend('RPE','NRPEstd','TPEmax','TPEmin','NRPEmax','NRPEmin','location','Best')
hold off

grid on;
xlabel('Sector No.');
ylabel(sprintf('PES (%s)', unit));
title(sprintf('NRPE6\\sigma rms=%.2f, RPEpp=%.2f, TPEpp=%.2f (N=%d)',...
	      Val.NRPE6sigma, Val.RPEpp, Val.TPEpp, rd));
set(gca, 'XLim', [0 Sct(end)])
%set(gca, 'YTick', [-200:20:200])
set(gca, 'XTick', [-100:20:300])
