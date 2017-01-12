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
% PN Code Parameters
N = 8;
PRBS = round(rand(1,N));
fc = 1.5e9;
Tp = N/fc;
% Carrier Parameters
c_fo = 79e9;
c_Fs = 2.1*c_fo;                  % Sampling frequency
c_T = 1/c_Fs;                     % Sampling period
c_L = c_Fs*Tp;                    % Length of signal in samples
c_t = (0:c_L-1)*c_T;              % Time vector
% define a scaling factor to expand array size of DSSS code to carrier size
sf = c_L/N;
% Define 
c = cos(2*pi*c_fo*c_t);           
c2 = cos(2*pi*c_fo*c_t+pi);
%% Scale DSSS array appropriately
for k = 1:N;
    if k == 1
        PRBS_new(1:1:round(sf)) = PRBS(1,1);
    else
        PRBS_new((k-1)*round(sf)+1:1:k*round(sf)) = PRBS(1,k);
    end
end

%% Mix the signals to create DSSS signal
dsss_sig = [];
for k=1:c_L
    if PRBS_new(1,k) == 0
        dsss_sig = [dsss_sig c];
    else
        dsss_sig = [dsss_sig c2];
    end
end
%%
dsss_n = 2^nextpow2(size(dsss_sig,2));
dsss_spectrum = fft(dsss_sig,dsss_n);
dsss_spectrum_mag = abs(dsss_spectrum/dsss_n);
dsss_spectrum_s = dsss_spectrum_mag(:,1:dsss_n/2+1);
dsss_spectrum_s(:,2:end-1) = 2*dsss_spectrum_s(:,2:end-1);
plot(0:(c_Fs/dsss_n):(c_Fs/2-c_Fs/dsss_n),dsss_spectrum_s(1:dsss_n/2))



