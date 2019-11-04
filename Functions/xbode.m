function [freq] = xbode(...
    sys_siso,...
    fig_color,...
    line_style,...
    fig_line_width,...
    freq_Hz,...
    linearAxis,...
    absMag)
% function xbode(sys_siso,fig_color,line_style,fig_line_width)
%
% =========================================================================
%   Copyright (c) 2008-, Xu Chen, xuchen@cal.berkeley.edu
%   Author(s): Xu Chen
% =========================================================================
% Initial Version: 2011-04-26
ni = nargin;
no = nargout;
error( nargchk(1, 7, ni) );
if nargin < 7
    absMag = 0;
end
if nargin <6
    linearAxis = 0;
end
if nargin <5 || isempty(freq_Hz)
    SW_FREQ = 0;
else
    SW_FREQ = 1;
end
if nargin <4
    fig_line_width = 1.5;
end
if nargin <3
    line_style = '-';
end
if nargin <2
    fig_color = 'b';
end

if fig_line_width < 0.5
    fig_line_width = 0.5;
end

[Mnum,Nnum,MMnum] = size(sys_siso);
fig_num = gcf;
for mm = 1:MMnum
    if mm == 1
        if SW_FREQ
            try
                [mag,ph,w]=bode(sys_siso(1,1,mm),freq_Hz*2*pi);
            catch
                [mag,ph,w]=bode(ss(sys_siso(1,1,mm)),freq_Hz*2*pi);
            end
        else
            try
                [mag,ph,w]=bode(sys_siso(1,1,mm));
            catch
                [mag,ph,w]=bode(ss(sys_siso(1,1,mm)));
            end
        end
    else
        try
            [mag,ph,w]=bode(sys_siso(1,1,mm),freq*2*pi);
        catch
            [mag,ph,w]=bode(ss(sys_siso(1,1,mm)),freq*2*pi);
        end
    end
    if mm == 1
        f = w/(2*pi);
        %         FreqMin=10^fix(log10(min(f))); % 2011-04-26
        FreqMin = min(f);
        FreqMax=max(f);
        freq = f; % for next loop
    end
    ph        = mod(ph+180,360)-180;
    magnitude = mag(1,1,:);
    magnitude = magnitude(:);
    phase     = ph(1,1,:);
    phase     = phase(:);
    if absMag
        mag2plot = magnitude;
    else
        mag2plot  = 20*log10(abs(magnitude));
    end
    %     fig_color = num2color(line_color_num, max([MMnum,max_line_color]));
    subplot(211);
    if linearAxis
        h1 = plot(f,mag2plot,...
            'LineStyle',line_style,...
            'Color',fig_color,'LineWidth',fig_line_width);
    else
        h1 = semilogx(f,mag2plot,...
            'LineStyle',line_style,...
            'Color',fig_color,'LineWidth',fig_line_width);
    end
    grid on;
    hold on;
    ylabel('Gain (dB)');
    %      xlabel('Frequency (Hz)');
    %      axis('off')
    subplot(212);
    if linearAxis
        h2 = plot(f,phase,...
            'LineStyle',line_style,...
            'Color',fig_color,'LineWidth',fig_line_width);
        
    else
        h2 = semilogx(f,phase,...
            'LineStyle',line_style,...
            'Color',fig_color,'LineWidth',fig_line_width);
    end
    grid on;
    hold on;
    ylabel('Phase (degree)');
    xlabel('Frequency (Hz)');
    %      line_color_num=line_color_num+1;
    %      axis([h1 h2],'equal')
end
figure(fig_num)
subplot(211);
grid on;
if absMag
    ylabel('Gain');
else
ylabel('Gain (dB)');
end
set( gca, 'xlim',  [FreqMin FreqMax]);
subplot(212);
grid on;
ylabel('Phase (degree)');
xlabel('Frequency (Hz)');
set( gca, 'ytick', [-180 -90  0 90 180] );
set( gca, 'ylim',  [-180 180]);
set( gca, 'xlim',  [FreqMin FreqMax]);
subplot(211)