function hd = ideal_lp(wc, M);
% M is length of the filter
% alpha gets the midpoint
% m is time-shifted version of 0..M-1
% wc   = radians (-pi..pi)

alpha = (M-1)/2; n = [0:1:(M-1)];
m = n - alpha;

fc = wc/pi;

% Hd(ejw) = {1.e-jaw  |w| <= wc
%           {0,       wc < |w| < pi

% hd      = sin(wc(n-a))/(pi*(n-a))
% sinc    = sin(pi*x)/(pi*x)

hd = fc*sinc(fc*m);