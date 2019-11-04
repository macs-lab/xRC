function TMR = getTMR(PES)
% 2010-01-18
% =========================================================================
%   Copyright 2008-, Xu Chen. All rights reserved
%   Author(s): Xu Chen, xuchen@cal.berkeley.edu
% =========================================================================

SIGMA   = sqrt(mean(PES.^2));
TMR     = 3*SIGMA;
