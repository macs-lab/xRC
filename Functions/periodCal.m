function [N,p] = periodCal(w,Ts,epsilon)
% see EE123 for periodic signal period calculation
%   Copyright (c) 2008-, Xu Chen
% Author(s): Xu Chen 



M = 1;
if w == 0
    N = 1;
    p = 0;
    else
    while (abs(floor(M*2*pi/w)-M*2*pi/w)>epsilon)
        M = M+1;
    end

    N = floor(M*2*pi/w);
    p = Ts*M*2*pi/w;
end
