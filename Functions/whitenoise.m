function [out] = whitenoise(N,Ts)
%WHITENOISE  Make a whitenoise with a unity power
%   function [out] = whitenoise(N,Ts)
%   Make a whitenoise with a unity power.
%   Input
%     N  - number of noise
%     Ts - sampling time
%   Output
%     out - noise
%
%   Copyright (c) 2004-2005, MSS benchmark working group
%   Copyright (c) 2006-, HDD benchmark working group
%   Ver.1.0, 2006-04-27

% Authors(s): M.Hirata
% Ver.1.0 2006-04-27

out = randn(N,1)/sqrt(2*Ts);

%% EOF of whitenoise.m
