%% Analog and Digital: Are we in sync about the sinc?
% Rect & Sinc: Analog and digital
%

%%
% The rect function is a very useful tool when analyzing signals 
% around us. It represents a signal that has a value of one for a given
% duration centered around time zero, and is zero everywhere else. 
N = 2^7;
M = 2^4-1;

ts = [0:1/N:1-1/N] - 1/2;
rect_vals = zeros(N,1);
rect_vals(N/2-(M-1)/2:N/2+(M+1)/2) = 1;
plot(ts, rect_vals, 'LineWidth', 2)
axis([-0.5, 0.5, -0.1, 1.1])

%% 
% We need a frequency domain understanding of the rect function 
% in order to fully understand it's behavior. In our early classes on signals 
% and systems, we were taught that the fourier transform of the rect 
% function is the sinc function. 
% 
% $$sinc(\omega \tau/2) = \tau \frac{\sin{(\omega \tau/2)}}{\omega \tau/2}$$
% 

%% 
% In practice when using digital computers, we are limited to 
% discrete-time representation of the rect function. 
% How then do we analyze this discrete-time signal in the frequency domain?

%%
% The DTFT or the discrete-time rect function provides a frequency 
% domain understanding of this signal, and mathematically, it can be shown 
% to be as follows:
%
% 
% $$X_{DTFT}(f) = \frac{\sin(\pi f M)}{\sin(\pi f)}.$$
% 

%% 
% This looks close enough to what MATLAB defines as sinc
% function which is given by:
%%
% 
% $$sinc(x) = \frac{\sin{\pi x}}{\pi x}$$
% 
%%
% Close, but not exactly the same! There is a pesky sin function in the
% denominator of the DTFT. Wouldn't it be nice if we could simply use the
% MATLAB sinc function and call it a day?

%%
% So, what are we giving up in making approximating the DTFT of the
% discrete-time rect by the sinc?

%% 
%%
% The fast Fourier transform or the FFT which is an efficient of the
% discrete Fouriet transform can be used for this purpose. 

% The DFT assumes that the signal is periodic. We use half pulses at
% beginning and end of the sequence to model a pulse that is symmetric
% around the zeroth sample. 

t = 0:1/N:(1-1/N);

x = zeros(N,1);
x(1:floor(M/2)+1) = 1; x((N-floor(M/2)+1):end) = 1;

%% Discrete-time Simulation
f = 0:1e-4:0.5;
% pi*f*M = 2*pi*(f/fs)*M/2  = w * tau/2
XDTFT_form = sin(pi*f*M)./sin(pi*f);

% f goes from [0, fs/2]
% sin(pi*f*M)/(pi*f) = (M * 1/fs)* sin(2*pi*(f/fs)*M/2)/(2*pi*(f/fs)*(M/2))
sinc_f = M*sinc(f*M);

figure;
plot(t, x, 'o-');

XDFT = fft(x);

%%
% The accuracy lost in this approximation is a result of aliasing caused by
% sampling the continuous rect function. The figure below compares
% the FFT of the discrete-time rect with formulas for the DTFT and pure
% sinc functions.


figure;
subplot(2, 1,1)
plot(1/N*[1:N/2], abs((XDFT(1:N/2))), '*-');
hold on
plot(f, abs(XDTFT_form), 'r')
plot(f, abs(sinc_f), 'm')
legend('FFT computed', 'DFT formula', 'Sinc formula')

subplot(2, 1,2);
plot(1/N*[1:N/2], abs((XDFT(1:N/2))), '*-');
hold on
plot(f, abs(XDTFT_form), 'r')
plot(f, abs(sinc_f), 'm')
axis([0.45, 0.5, -inf, 2])
grid on;
legend('FFT computed', 'DFT formula', 'Sinc formula')

% figure;
% plot(1/N*[1:N/2], unwrap(angle((XDFT(1:N/2)))));

% figure;
% subplot(2,1,1);
% plot(1/N*[1:N/2], real(XDFT(1:N/2)), '*-')
% subplot(2,1,2);
% plot(1/N*[1:N/2], imag(XDFT(1:N/2)),'*-');
% axis([0, 0.5, -pi, pi])

%%
% At normalized frequency of 0.5, the approximation
% error could be about 36%.

%% Sinc aliasing
falias = [-40.5:1:40.5];
sincout = M*sinc(-falias*M);
sincalias = sum(sincout)
disp(['Sinc Formula = ', num2str(sinc_f(end)), '  DTFT = ', num2str(sincalias)])
approx_error = ((XDTFT_form(end)-sinc_f(end))/XDTFT_form(end))*100
disp(['Approximation error = ', num2str(approx_error), '%'])

%% 
% A way to visualize aliasing that is occurring is by drawing
% frequency-shifted but overlapping sinc functions. The 
% response at a particular frequency is the sum of these sinc's.
figure; hold on;
f = -100:1e-3:100;
sinc_at_fcenter = M*sinc((f)*M);
fcenter = [-40:40];
for ii=1:length(fcenter)
    plot((f+fcenter(ii)), sinc_at_fcenter, 'LineWidth', 1.5)
    axis([-1, 1, -5, M+3])
end

%% 
% So, analog or digital, the sinc function is incredibly useful.
% Just keep in mind that the aliasing in the discrete-time world causes it
% to be somewhat of an approximation. Aliasing, which results in a number 
% frequency-shifted overlapping sinc functions, creates a $\sin(\pi x)$
% term in the denominator of the discrete-time Fourier transform, in place
% of the $\pi x$ in the sinc function.

%%
% Analog or digital, left-handed or right-handed, we can all continue to
% embrace the rect function to better understand the world around
% us. 

%%
% Kumbaya!
    
