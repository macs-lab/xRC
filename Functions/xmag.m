function varargout = xmag(...
    sys_siso,...
    fig_color,...
    line_style,...
    fig_line_width,...
    freq_Hz,...
    sw_absMag)
% function xbode(sys_siso,fig_color,line_style,fig_line_width)
%
% Author: Xu Chen xuchen@cal.berkeley.edu
% Initial Version: 2011-04-26
ni = nargin;
no = nargout;
error( nargchk(1, 6, ni) );
if ni < 6
    sw_absMag = 0;
end
if nargin <5
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
    ph          = mod(ph+180,360)-180;
    magnitude   = mag(1,1,:);
    magnitude   = magnitude(:);
    phase       = ph(1,1,:);
    phase       = phase(:);
    
    if sw_absMag
        h1 = plot(f,abs(magnitude),...
            'LineStyle',line_style,...
            'Color',fig_color,'LineWidth',fig_line_width);    
    else
        h1 = semilogx(f,20*log10(abs(magnitude)),...
            'LineStyle',line_style,...
            'Color',fig_color,'LineWidth',fig_line_width);
    end
    grid on;
    hold on;
    if sw_absMag
        ylabel('Magnitude (abs)');
        set(gca, 'ylim', [0.8*min(magnitude),max(magnitude)*1.5])
    else
        ylabel ('Magnitude (dB)')
        minaxis = min(20*log10(magnitude));
        maxaxis = max(20*log10(magnitude));
        if minaxis == 0
            minaxis = -5;
        end
        if maxaxis == 0
            maxaxis = 5;
        end
%         if abs((maxaxis - minaxis)/maxaxis)> || 
        if (maxaxis - minaxis)>20
            set(gca, 'ylim', [0.8*minaxis,maxaxis*1.5])
        end
    end
    xlabel('Frequency (Hz)');
    %      line_color_num=line_color_num+1;
    %      axis([h1 h2],'equal')
end
set( gca, 'xlim',  [FreqMin FreqMax]);

if nargout == 1
    varargout{1} = freq;
elseif nargout == 2
    varargout{1} = freq;
    if sw_absMag
        varargout{2} = magnitude;
    else
        varargout{2} = 20*log10(magnitude);
    end
end