% Xu Chen xuchen@cal.berkeley.edu  2010-01-18
%% baseline disturbance define
% load mainDistDataV1;
% revo = 30;
% ============================================================
%   Copyright (c) 2008-, Xu Chen, chx@uw.edu
%   Author(s): Xu Chen
% ============================================================
% 2010-11-17 :
%           added arbitrary user defined number of narrow band
%           components
folder = fileparts(which(mfilename));
addpath(genpath(folder));

disp('=====================================================')
disp(' Input revo: (20 for default)')
revo = input('');
if isempty(revo)
    revo            = 20;
end
Tsim            = revo*60/PlantData.rpm; %   simulation time.
[FlutterDist,SensorNoise, ForceDist] = SetDisturbance(Tsim,PlantData);
RRO             = SetRRO(Tsim,PlantData);
% DesignMethod & DesignerName
strDesignMethod = 'PID + Notch Filter + DOB + Adaptive Estimation';
strDesigner     = 'Xu Chen';
% Excel_filename =[pwd,'\eval_sheet.xls'];

Ts              = PlantData.Ts; %   sampling time of PES
% Tc = Ts/2; % sampling time of control input (Multi-rate Notch Filter)
Tc              = Ts;

if ~exist('SW_VARIABLE_DIST_BAND','var') || SW_VARIABLE_DIST_BAND == 0
    % \\\\\\\\\\\\\\\narrow band disturbances define
    if ~exist('NBn','var')
        NBn             = 2;%               number of narrow bands
    end
    % NBn = 3;
    % NBn = 4;
    
    if NBn==2
        % NBw = [0.5/6*pi/Ts 0.25/6*pi/Ts];
        %                               narrow band frequencies
        
        % NBw = [400*2*pi 650*2*pi];
        %                               close frequency
        
        % NBw = [500*2*pi 1000*2*pi];
        %                               harmonic case
        
        NBw = [500*2*pi 1200*2*pi];
        %                               1100 very noisy there
        
        % NBw = [500*2*pi 0*2*pi];
        %                               one component but n = 2 in adaptation
        
%         NBw = [800*2*pi 1100*2*pi];
        %                               one freq higher than the BW
        
        % NBw = [504*2*pi 648*2*pi];
        %                               experiment data from Guoxiao's paper
        
        % NBw = [504*2*pi 696*2*pi];
        %                               experiment data from Guoxiao's paper
        
        NBw2 = NBw;
        
        % NBw2 = [350*2*pi 900*2*pi];
        %                               one freq higher than the BW
        
        % NBw = NBw2;
        % NBw = [500*2*pi 1800*2*pi];
        %                               one freq higher than the BW
        
        NBamp = [4e-4 3e-4];
        %                               dist amplitude
        
        % NBamp = [2e-4 1.5e-4];
        %                               dist amplitude
    elseif NBn==3
        NBw = [504*2*pi 800*2*pi 1200*2*pi];
        %                               three disturbances
        % NBw = [504*2*pi 648*2*pi 1200*2*pi];
        %                               three disturbances
        NBamp = [4e-4 3e-4 4e-4];
        %                               dist amplitude
    elseif NBn==4
        NBw = [504 650 800 900]*2*pi;
        %                               four disturbances
        NBamp = [4e-4 3e-4 4e-4 3e-4];
        %                               dist amplitude
        % NBw = [504*2*pi 648*2*pi 1200*2*pi];
        %                               three disturbances
    elseif NBn == 1
        NBw         = 1000*2*pi;
        NBamp       = 7e-4;
    elseif NBn == 5
        NBw         = [504 650 800 900 1200]*2*pi;
        NBamp       = ones(NBn,1)/NBn*7e-4;
    end
    try
        f1          = NBw(1)/2/pi;
        f2          = NBw(2)/2/pi;
        % NBamp = [1 0.3];
        %                                   dist amplitude
        [w1N, w1period]     = periodCal(NBw(1)*Ts,Ts,0.001);
        [w2N, w2period]     = periodCal(NBw(2)*Ts,Ts,0.001);
    catch
    end
elseif SW_VARIABLE_DIST_BAND == 1
    % \\                                load the last seed from rand
    %                                   ensuring different rand results
    RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*cputime)));
    iter_temp = 0;
    while(1)
        iter_temp = iter_temp +1;
        if exist('NBf_upperBound','var') && NBf_upperBound ~= 0
            if ~exist('NBf_lowerBound','var') 
                NBf_lowerBound = 300;
            end
            NBw         = 2*pi* ...
                ( rand(1,NBn) * (NBf_upperBound-NBf_lowerBound) ) ...
                + 2*pi* NBf_lowerBound;
            if NBn >= 3
                NBf_gap     = 150;
            else
                NBf_gap     = 180;
            end
        else
            NBw         = 2*pi* (rand(1,NBn)*900) + 2*pi* 300;
            NBf_gap     = 250;
        end
        
        if NBn > 1
            if min(abs(diff(NBw/2/pi))) >= NBf_gap
                break;
            else
                if iter_temp > 50
                    error('The disturbance frequency region may be too small')
                end
            end
            if NBn == 2
                [w1N, w1period]     = periodCal(NBw(1)*Ts,Ts,0.001);% epsilon = 0.001
                [w2N, w2period]     = periodCal(NBw(2)*Ts,Ts,0.001);
            end
        else % NBn == 1
            break;
        end
    end
    NBamp       = ones(1,NBn)/NBn*7e-4; 
    %                                   dist amplitude
elseif SW_VARIABLE_DIST_BAND == 2
    % manual frequency input
    if ~exist('NBn','var')
        NBn             = 2;%               number of narrow bands
    end    
    commandwindow
    for mm = 1:NBn
        disp([' Input the ',num2str(mm),'-th ', 'frequency in Hz:'])
        NBf(mm) = input('');
    end
    NBw = NBf*2*pi;
    if NBn==2
        NBw2 = NBw;
        NBamp = [4e-4 3e-4];
        %                               dist amplitude
    elseif NBn==3
        NBamp = [4e-4 3e-4 4e-4];
        %                               dist amplitude
    elseif NBn==4
        NBamp = [4e-4 3e-4 4e-4 3e-4];
        %                               dist amplitude
    elseif NBn == 1
        NBamp       = 7e-4;
    elseif NBn == 5
        NBamp       = ones(NBn,1)/NBn*7e-4;
    end
    try
        f1          = NBw(1)/2/pi;
        f2          = NBw(2)/2/pi;
        % NBamp = [1 0.3];
        %                                   dist amplitude
        [w1N, w1period]     = periodCal(NBw(1)*Ts,Ts,0.001);
        [w2N, w2period]     = periodCal(NBw(2)*Ts,Ts,0.001);
    catch
    end
else
end

NBm             = 0;
%                                   first m revolution no NB dist
t_NBon          = NBm/120;
%                                   NB Dist injection time
t_Qon           = (NBm+3-NBm)/120;
%                                   Peak Q filter on time
% t_AdapOn = (NBm+1)/120;
%                                   Adaptation start time
% t_AdapOn = t_NBon;
t_AdapOn        = 0/120;
% t_Qon = Ts;
NBf             = NBw/2/pi;
%                                   freq in Hz
NBphi           = rand(length(NBw),1);
%                                       dist phase
MNDist.time     = [0:(revo*PlantData.num_servo-1)]'*Ts; 
%                               2011-04-27 added "-1"

if isempty(flag_NB_method)
    MNDist.data     = signalgen_ti(...
        Ts,NBw,NBamp,NBphi,NBn,revo*PlantData.num_servo-1,NBm,...
        'noplot','noiseOff');
    %                               2011-04-27 added "-1"
    %     MNDist.data = signalgen_tvf(...
    %         Ts,NBw,NBw2,NBamp,NBphi,NBn,revo*PlantData.num_servo,NBm,'plot','noiseOff');
else
%     MNDist.data     =
%     NBgen_white_excite(NBw,revo*PlantData.num_servo+1,Ts);
    MNDist.data     = NBgen_white_excite(NBw,revo*PlantData.num_servo,Ts);
%     2011-04-27 

end
MNDist.Ts           = Ts;
z                   = tf('z',Ts);
% 
% specDist = specCal(MNDist.data,1/Ts);
% figure,plot(specDist.f,specDist.amp,bode_opt)
% title 'dist spec'
