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
f = fs*linspace(0,1,NFFT/2+1);
P2 = abs(outfft/NFFT);
P1 = P2(:,1:NFFT/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1);
plot(f,P1)
% plot(f,y_mag(1:size(f,2))/size(f,2));
% title('Frequency Response')
% xlabel('Frequency [Hz]');


%%
clear all
close all
clc
N = 1024;
PRBS = round(rand(1,N));

c_fo = 79e9;
c_Fs = 3*c_fo;                    % Sampling frequency
c_T = 1/c_Fs;                     % Sampling period
c_L = 161792;                     % Length of signal in samples
c_t = (0:c_L-1)*c_T;              % Time vector

c = cos(2*pi*c_fo*c_t);          % First row wave
c2 = cos(2*pi*c_fo*c_t+pi);
% c_n = 2^nextpow2(c_L);
% c_spectrum = fft(c,c_n);
% c_spectrum_mag = abs(c_spectrum/c_n);
% c_spectrum_s = c_spectrum_mag(:,1:c_n/2+1);
% c_spectrum_s(:,2:end-1) = 2*c_spectrum_s(:,2:end-1);
% plot(0:(c_Fs/c_n):(c_Fs/2-c_Fs/c_n),c_spectrum_s(1:c_n/2))

% Scale PBRS array appropriately
% each sample bit gets expanded into 158 bits
%PRBS_new(1:1:159) = PRBS(1,1);
%PRBS_new(160:1:318) = PRBS(1,2);
%PRBS_new(319:1:477) = PRBS(1,3); etc

for k = 1:N;
    if k == 1
        PRBS_new(1:1:158) = PRBS(1,1);
    else
        PRBS_new((k-1)*158+1:1:k*158) = PRBS(1,k);
    end
end

%%
dsss_sig = [];

for k=1:c_L
    if PRBS(1,k) == 0
        dsss_sig = [dsss_sig c];
    else
        dsss_sig = [dsss_sig c2];
    end
end
%%

dsss_n = 2^nextpow2(size(dsss_sig));
dsss_spectrum = fft(dsss_sig,dsss_n);
dsss_spectrum_mag = abs(dsss_spectrum/dsss_n);
dsss_spectrum_s = dsss_spectrum_mag(:,1:dsss_n/2+1);
dsss_spectrum_s(:,2:end-1) = 2*dsss_spectrum_s(:,2:end-1);
plot(0:(c_Fs/dsss_n):(c_Fs/2-c_Fs/dsss_n),dsss_spectrum_s(1:dsss_n/2))

% dsss_size = size(dsss_sig, 2);
% NFFT = 2^nextpow2(dsss_size); % Next power of 2 from length of y
% DSSS_spectrum = fft(dsss_sig,NFFT);
% f = 81.5e9*linspace(0,1,NFFT/2+1);
% DSSS_P2 = abs(DSSS_spectrum/NFFT);
% DSSS_P1 = DSSS_P2(:,1:NFFT/2+1);
% DSSS_P1(:,2:end-1) = 2*DSSS_P1(:,2:end-1);
% plot(f,DSSS_P1)


%figure,plot([1:12120],dsss_sig);
%axis([-1 12220 -1.5 1.5]);
%title('\bf\it DSSS Signal');

% Plotting the FFT of DSSS signal
%figure,plot([1:12120],abs(fft(dsss_sig)))