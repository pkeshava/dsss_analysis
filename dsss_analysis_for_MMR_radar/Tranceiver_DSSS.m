% Radar Tranceiver Spread Spectrum
clear all 
close all
clc

% Initialize parameters
% Define BW/2 (fc) a value that wil give a achievable range resolution of
% 0.1m
fc = 1.5e9;
fs = 2*fc;          % Since signal highest frequency element is fc then we sample at 2*fc
tc = 1/fc;
% Define the number of bits in the code to give a suported unambigous range
% of Rmax ~100m
N = 1024;

Tp = N*tc;                  % Length of signal in seconds
L = Tp*fc;                  % Length of signal in samples
t = (0:L-1)*tc;             % Time vector


PRBS = round(rand(1,N));
plot(t,PRBS)
axis([0 7e-7 -.5 1.5]);
title('Pseudorandom Spreading Code For Transceiver')
xlabel('Time [s]');
ylabel('Value');


%%
% Define Y as FFT of X
% Y = fft(PRBS);
M = size(PRBS, 2);
NFFT = 2^nextpow2(M); % Next power of 2 from length of y
outfft = fft(PRBS,NFFT);
f = fs/2*linspace(0,1,NFFT/2+1);
P2 = abs(outfft/NFFT);
P1 = P2(:,1:NFFT/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1);
plot(f,P1)
% plot(f,y_mag(1:size(f,2))/size(f,2));
% title('Frequency Response')
% xlabel('Frequency [Hz]');


%%
dsss_sig = [];
   
fo = 79e9;
c1 = cos(2*pi*fo*t);
c2 = cos(2*pi*fo*t+pi);

for k=1:N
    if PRBS(1,k) == 0
        dsss_sig = [dsss_sig c1];
    else
        dsss_sig = [dsss_sig c2];
    end
end

M = size(dsss_sig, 2);
NFFT = 2^nextpow2(M); % Next power of 2 from length of y
DSSS_spectrum = fft(dsss_sig,NFFT);
f = 158e9/2*linspace(0,1,NFFT/2+1);
DSSS_P2 = abs(DSSS_spectrum/NFFT);
DSSS_P1 = DSSS_P2(:,1:NFFT/2+1);
DSSS_P1(:,2:end-1) = 2*DSSS_P1(:,2:end-1);
plot(f,DSSS_P1)


%figure,plot([1:12120],dsss_sig);
%axis([-1 12220 -1.5 1.5]);
%title('\bf\it DSSS Signal');

% Plotting the FFT of DSSS signal
%figure,plot([1:12120],abs(fft(dsss_sig)))