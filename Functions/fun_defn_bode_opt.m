function bode_opt = fun_defn_bode_opt(xlim_low_Hz, Ts, flag_linear)
% function bode_opt = fun_defn_bode_opt(xlim_low_Hz, Ts, flag_linear)
% =========================================================================
%   Copyright(C) 2008-, Xu Chen, xuchen@cal.berkeley.edu
%   Author(s): Xu Chen
%   All rights reserved
% =========================================================================
% first ver: 2011-01-03
if nargin < 3
    flag_linear = 'logScale';
    if nargin < 2
        Ts = 1/26400;
        if nargin < 1
            xlim_low_Hz = 100;
        end
    end
end
if strcmp(flag_linear,'logScale') || strcmp(flag_linear,'log_scale')
    bode_opt = bodeoptions;
    try
        bode_opt.Xlim = [xlim_low_Hz, 0.5/Ts];
    catch
        bode_opt.Xlim{1} = [xlim_low_Hz, 0.5/Ts];
    end
    bode_opt.FreqUnits = 'Hz';
    bode_opt.PhaseWrapping = 'On';
else
    bode_opt = bodeoptions;
    try
        bode_opt.Xlim = [xlim_low_Hz, 0.5/Ts];
    catch
        bode_opt.Xlim{1} = [xlim_low_Hz, 0.5/Ts];
    end
    bode_opt.FreqUnits = 'Hz';
    bode_opt.FreqScale = 'linear';
    bode_opt.MagScale = 'linear';
    bode_opt.PhaseWrapping = 'On';
end