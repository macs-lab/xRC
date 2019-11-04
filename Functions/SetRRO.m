function [RRO] = SetRRO(t_simu, PlantData);
%SetRRO  Set repeatable runout (RRO)
%   [RRO] = SetRRO(PlantData,DistParam);
% ============================================================
%   Copyright (c) 2008-, Xu Chen, xuchen@cal.berkeley.edu
%   Author(s): Xu Chen 
% ============================================================



Ts = PlantData.Ts;
Ty = PlantData.Ty;
Tu = PlantData.Tu;

rpm = PlantData.rpm; %5400;
num_servo = PlantData.num_servo;

% Data from plant model
Fns = 1/Ts/2; % Nyquist freq. of PES sampling
Fnu = 1/Tu/2; % Nyquist freq. of plant input
Fny = 1/Ty/2; % Nyquist freq. of plant output

dFs = 6; % Frequency resolution of PES
dFu = 6; % Frequency resolution of plant input
dFy = 6; % Frequency resolution of plant output
Nums  = Fns/dFs;
% Numu  = Fnu/dFu;
% Numy  = Fny/dFy;
% Freqs = dFs*(1:Nums)';
% Frequ = dFu*(1:Numu)';
% Freqy = dFy*(1:Numy)';
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

FreqRRO     = [ 1,  2,  3]*(rpm/60);  % Frequency [Hz]
AmpRRO      = [0.15, 0.03, 0.006];    % unit:[Track] at PES, Amplitude of each RRO component


% seeds for random signal
Seed_RRODist      = 1;
% Frequency vector is defined by PES sampling frequency
Freq = dFs*(1:Nums)';
NUM  = Fns/dFs;

% RRO dist. in time domain
NUM_RRODist = t_simu/Ty+1;
% Check RROSequence
load RROSequence
RROSequence_rep = ...
    [repmat(RROSequence(1:PlantData.num_servo),t_simu/Ts/num_servo,1);0];
% different from the original HDDBenchmark since this is not track following

rand('state',Seed_RRODist); % Here is rand command (not randn)
[num_RRODist] = length(FreqRRO);
RRODist  = zeros(NUM_RRODist,1);
RRODistT   = (0:NUM_RRODist-1)'*Ts;

for ii=1:num_RRODist
    Frro = FreqRRO(ii);
    Arro = AmpRRO(ii);
    ph_lag = rand(1)*2*pi;
    RRODist = RRODist + sin(2*pi*Frro*RRODistT + ph_lag)*Arro; % fixed (Ver.1.1)
end
RRODist = RRODist + RROSequence_rep;

% RRO dist. in frequency domain
Spec    = zeros(NUM,1);
[N_RRO] = length(FreqRRO);
cc      = 0;
for(jj=1:N_RRO)
    ii = find(Freq == FreqRRO(jj));
    if(ii > 0)
        Spec(ii(1)) = AmpRRO(jj);
        cc = cc + 1;
    end
end
if ( cc ~= N_RRO )
    error('Frequency vector Freqs is not correct')
end

% Output parameters in time domain
RRO.Data = RRODist;
RRO.Time = RRODistT;
RRO.Ts   = Ts;
% Output parameters in frequency domain
RRO.Freq = Freq;
RRO.Spec = Spec;

%% EOF of SetRRO.m
