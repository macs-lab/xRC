function varargout = lattice_filter(fc, bw, Ts)
% Design lattice based notch filter G
%   G = (1+A)/2
% where A is a digital all-pass filter, enjoying the fewest number of
% multiplication
% Usage: sys_notch = notch_lattice(fc, bw)
% where fc = center frequency of notch in Hz
%       bw = 3-dB bandwidth of the notch filter
%       Ts = sampling time
% uses LTI object provided by LTI toolbox
% Date: 7/17/08
% Author: Jake Y-H Kim

w0 = 2*pi*fc*Ts;
Omega = 2*pi*bw*Ts;
k1 = -cos(w0);
k2 = (1-tan(Omega/2))/(1+tan(Omega/2));
if 0
    z = tf('z',Ts);
    sys_A = (k2 + k1*(1+k2)*z^-1 + z^-2)/(1 + k1*(1+k2)*z^-1 + k2*z^-2);
else
    th2 = asin(k2);
    th = w0;
    A = [cos(th) -sin(th2)*sin(th)
        sin(th)  sin(th2)*cos(th)];
    B = [cos(th2)*sin(th); -cos(th2)*cos(th)];
    C = [0 cos(th2)];
    D = sin(th2);
    sys_A = ss(A,B,C,D,Ts);
end

G = minreal((1+sys_A)/2);       % notch
H = minreal((1-sys_A)/2);       % complementary BPF
CosSinFw1 = round(2^15*[cos(th) sin(th)]);
CosSinFw2 = round(2^15*[cos(th2) sin(th2)]);
FW.CosSin1 = CosSinFw1;
FW.CosSin2 = CosSinFw2;

if nargout==0
    bode(G,H)
    fprintf('k1 = %4.2f, k2 = %4.2f\n', k1, k2);
    fprintf('(cos(th1), sin(th1))_FW = (%5i, %5i)\n', CosSinFw1(1), CosSinFw1(2));
    fprintf('(cos(th2), sin(th2))_FW = (%5i, %5i)\n', CosSinFw2(1), CosSinFw2(2));
else
    varargout{1} = G;
    varargout{2} = H;
    varargout{3} = FW;
end
