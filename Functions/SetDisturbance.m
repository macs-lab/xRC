function [FlutterDist, SensorNoise, ForceDist] = SetDisturbance(t_simu,PlantData);
% Disturbance and noise profile for spiral writing
%% disturbance parameters
% ============================================================
%   Copyright (c) 2008-, Xu Chen, xuchen@cal.berkeley.edu
%   Author(s): Xu Chen 
% ============================================================



Ts = PlantData.Ts;
Ty = PlantData.Ty;
% Tu = PlantData.Tu;
Tu = Ts;

% Data from plant model
Fns = 1/Ts/2; % Nyquist freq. of PES sampling
Fnu = 1/Tu/2; % Nyquist freq. of plant input
Fny = 1/Ty/2; % Nyquist freq. of plant output

dFs = 6; % Frequency resolution of PES
dFu = 6; % Frequency resolution of plant input
dFy = 6; % Frequency resolution of plant output

% Check frequency resolution
if ( (Fns/dFs) ~= fix(Fns/dFs) )
    error('Frequency resolution dF0s is not appropriate')
end
if ( (Fnu/dFu) ~= fix(Fnu/dFu) )
    error('Frequency resolution dF0u is not appropriate')
end
if ( (Fny/dFy) ~= fix(Fny/dFy) )
    error('Frequency resolution dF0y is not appropriate')
end

Nums  = Fns/dFs;
Numu  = Fnu/dFu;
Numy  = Fny/dFy;
Freqs = dFs*(1:Nums)';
Frequ = dFu*(1:Numu)';
Freqy = dFy*(1:Numy)';

% Scaling parameter
scl = sqrt(2*Ts);

% Sensor noise
AmpSensorNoise = 1.5e-2*scl; 
%       unit:[Track], 1 sigma of Sensor Noise

% Force disturbance
AmpForceDist = 1.0e-4*scl; 
%       unit:[A] at control input, 1 sigma of Torque Noise with sampling Ts

FlutterFreq    = [  750,  780,  900, 1020, 1080, 1230, 1800, 3000, 5000];
% FlutterFreqSigma = [   40,   30,   30,   30,   40,   40,   40,   50,   50];
FlutterFreqZeta = [ 0.01, 0.01, 0.01, 0.01,0.005, 0.01,0.005,0.002,0.002];
AmpFlutterDist = [ 0.09, 0.17, 0.20, 0.17, 0.06, 0.09, 0.06, 0.06, 0.12]*scl;

% Seed for random signal, the same seed results in the same output from the
% rand function. if we don't set the seed, all dist will be the same.
Seed_ForceDist    = 2; 
Seed_SensorNoise  = 3;
Seed_FlutterDist  = 4;

% Frequency vector
DistParam.Nums              = Nums;
DistParam.Numu              = Numu;
DistParam.Numy              = Numy;
DistParam.dFs               = dFs;
DistParam.dFu               = dFu;
DistParam.dFy               = dFy;
DistParam.Freqs             = Freqs;
DistParam.Frequ             = Frequ;
DistParam.Freqy             = Freqy;
% Sensor noise
DistParam.AmpSensorNoise    = AmpSensorNoise;
% Force disturbance
DistParam.AmpForceDist      = AmpForceDist;
% Flutter disturbance
DistParam.AmpFlutterDist    = AmpFlutterDist;
DistParam.FlutterFreq       = FlutterFreq;
DistParam.FlutterFreqZeta   = FlutterFreqZeta;
% Seeds for random signal
DistParam.Seed_ForceDist    = Seed_ForceDist;
DistParam.Seed_SensorNoise  = Seed_SensorNoise;
DistParam.Seed_FlutterDist  = Seed_FlutterDist;
DistParam.Tu                = Tu;
DistParam.Ts                = Ts;
DistParam.Ty                = Ty;
DistParam.t_simu            = t_simu;

%% Data output
ForceDist = SetForceDist(DistParam);
SensorNoise = SetSensorNoise(DistParam);
FlutterDist = SetFlutterDist(DistParam);
%% Plot the disturbance in time and freq domain
% time domain
% figure;
% subplot(311)
% plot(FlutterDist.Time,FlutterDist.Data);
% title('Flutter Disturbance');ylabel('Track');
% subplot(312)
% plot(ForceDist.Time,ForceDist.Data);
% title('Force Disturbance');ylabel('Track');
% subplot(313)
% plot(SensorNoise.Time,SensorNoise.Data)
% title('Sensor Noise');xlabel('Time[sec]');ylabel('Ampere');
% 
% % frequency domain
% NRROSpec = sqrt(...
%   ForceDist.SpecAtPes.^2 + FlutterDist.Spec.^2 + SensorNoise.Spec.^2);
% SpecAll = NRROSpec;
% figure;
% % semilogx(...
% % 	Freqs,...
% %      20*log10(...
% %      [FlutterDist.Spec, ForceDist.Spec, SensorNoise.Spec, SpecAll]));
% semilogx(...
%   Freqs,...
%   20*log10(...
%   [FlutterDist.Spec, ForceDist.Spec, SensorNoise.Spec, SpecAll]));
% title('Dist data in frequency domain')
% xlabel('Frequency [Hz]'), ylabel('Magnitude [dB]');grid on;
% % axis([1e1 1e4 -120 0])
% legend('Fluttter Dist','Force Dist','Sensor Noise','ALL')

%% set sensor noise
    function [SensorNoise] = SetSensorNoise(DistParam);
        %SetSensorNoise  Set sensor noise in time and frequency domain

        % simulation parameters
        Ts                 = DistParam.Ts; % PES Sampling
        % seeds for random signal
        Seed_SensorNoise  = DistParam.Seed_SensorNoise;

        % Frequency vector is defined by PES sampling frequency
        Freq = DistParam.Freqs;
        NUM  = DistParam.Nums;

        % Sensor noise in time domain
        NUM_SensorNoise = DistParam.t_simu/Ty+1;
        randn('state',Seed_SensorNoise);
        Data = whitenoise(NUM_SensorNoise,Ts)*DistParam.AmpSensorNoise;
        Time = (0:NUM_SensorNoise-1)'*Ts;

        % Sensor noise in frequency domain
        Spec = DistParam.AmpSensorNoise*ones(NUM,1);

        % Output parameters in time domain
        SensorNoise.Data = Data;
        SensorNoise.Time = Time;
        SensorNoise.Ts   = Ts;
        % Output parameters in frequency domain
        SensorNoise.Freq = Freq;
        SensorNoise.Spec = Spec;
    end
%% EOF of SetSensorNoise.m

%% set force dist
    function [ForceDist] = SetForceDist(DistParam);
        % Read parameters
        AmpForceDist       = DistParam.AmpForceDist;
        Ts                 = DistParam.Ts; % PES Sampling
        Tu                 = DistParam.Tu; % Plant input
        % Seeds for random signal
        Seed_ForceDist  = DistParam.Seed_ForceDist;
        % Frequency vector is defined by PES sampling frequency
        Freq = DistParam.Freqs;
        NUM  = DistParam.Nums;

        % Force dist. in time domain
        if fix(Ts/Tu) ~= (Ts/Tu)
            error('Ts/Tu must be integer');
        end
        NUM_ForceDist = DistParam.t_simu/Tu+1;
        randn('state',Seed_ForceDist);
        Data       = whitenoise(NUM_ForceDist,Tu)*DistParam.AmpForceDist;
        Time       = (0:NUM_ForceDist-1)'*Tu;
        % dist at PES
        DistOut    = lsim(PlantData.Pn,Data,Time);
        dcidx      = 1:(Ts/Tu):NUM_ForceDist;
        DataAtPes  = DistOut(dcidx); % down sampling by (Ts/Tu)
        TimeAtPes  = (0:length(DataAtPes)-1)'*Ts;

        % Force dist. in frequency domain
        Spec      = DistParam.AmpForceDist*ones(NUM,1);
        [mag,phs] = bode( PlantData.Pn*DistParam.AmpForceDist, 2*pi*Freq );
        SpecAtPes = mag(:);

        % Output parameters in time domain
        ForceDist.Data      = Data;
        ForceDist.Time      = Time;
        ForceDist.Ts        = Tu;
        ForceDist.DataAtPes = DataAtPes;
        ForceDist.TimeAtPes = TimeAtPes;
        ForceDist.TsAtPes   = Ts;
        % Output parameters in frequency domain
        ForceDist.Freq      = Freq;
        ForceDist.Spec      = Spec;
        ForceDist.SpecAtPes = SpecAtPes;
    end
%% EOF SetForceDist.m

%% set flutter dist
    function [FlutterDist] = SetFlutterDist(DistParam);
        % SetFlutterDist  
%         Set flutter disturbance in time and frequency domain

        % Read parameters
        AmpFlutterDist     = DistParam.AmpFlutterDist;
        FlutterFreq        = DistParam.FlutterFreq;
        FlutterFreqZeta    = DistParam.FlutterFreqZeta;
        Ts                 = DistParam.Ts; % PES Sampling
        Ty                 = DistParam.Ty; % Plant output
        % Seeds for random signal
        Seed_FlutterDist  = DistParam.Seed_FlutterDist;
        % Frequency vector is defined by PES sampling frequency
        Freq = DistParam.Freqs;
        NUM  = DistParam.Nums;

        % sys_FlutterDist
        lg = length(AmpFlutterDist);
        for ii=1:lg
            w_f    = 2*pi*FlutterFreq(ii);
            zeta_f = FlutterFreqZeta(ii);
            pk_f   = AmpFlutterDist(ii);
            sys_FlutterDist(ii) = tf(...
                pk_f*2*zeta_f*w_f^2,[1 2*zeta_f*w_f w_f^2] );
        end

        % Flutter dist. in time domain
        if fix(Ts/Ty) ~= (Ts/Ty)
            error('Ts/Ty must be integer');
        end
        NUM_FlutterDist = DistParam.t_simu/Ty+1;
        randn('state',Seed_FlutterDist);
        Data = zeros(NUM_FlutterDist,1);
        Time = (0:NUM_FlutterDist-1)'*Ty; % Ver.2.0
        for ii=1:lg
            sys_FlutterDist0 = sys_FlutterDist(ii);
            FlutterDistIn = whitenoise(NUM_FlutterDist,Ty); % Ver.2.0
            Data = Data + lsim(sys_FlutterDist0,FlutterDistIn,Time);
        end
        dcidx     = 1:(Ts/Ty):NUM_FlutterDist;
        DataAtPes = Data(dcidx); % down sampling by (Ts/Ty)
        TimeAtPes = (0:length(DataAtPes)-1)'*Ts;

        % Flutter dist in frequency domain
        Spec2 = zeros(NUM,1);
        for ii=1:lg
            [mag,phs] = bode(sys_FlutterDist(ii), 2*pi*Freq );
            Spec2 = Spec2 + mag(:).^2;
        end
        Spec  = sqrt( Spec2 );

        % Output parameters in time domain
        FlutterDist.Data      = Data;
        FlutterDist.Time      = Time;
        FlutterDist.Ts        = Ty;
        FlutterDist.DataAtPes = DataAtPes;
        FlutterDist.TimeAtPes = TimeAtPes;
        FlutterDist.TsAtPes   = Ts;
        % Output parameters in frequency domain
        FlutterDist.Freq      = Freq;
        FlutterDist.Spec      = Spec;

    end
%% EOF of SetFlutterDist.m

end
%% EOF SetDisturbance.m
