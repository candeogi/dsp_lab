clearvars
close all
clc

%% LAB 6 DOWNSAMPLING

%% load audio

% taken from https://helpguide.sony.net/high-res/sample1/v1/en/index.html
[x,Fs] = audioread('Sample_BeeMoved_48kHz16bit.m4a'); % sampled at Fs = 48kHz
T = 1/Fs;
Nx = length(x);
F = 1/(Nx*T);

% listen to first 10 seconds
%sound(x(1:10*Fs,:),Fs)
X = T*fft(x);               
% show the absolute value of x in frequency
time = T*(0:Nx-1);
frequency = F*(0:Nx-1); 

%plot signal in time domain
figure(1);
plot(frequency/1e3,20*log10(abs(X))); grid;
xlabel('frequency [kHz]'); xlim([0 Fs/1e3]);
title('original audio signal in frequency');

%% filter audio

% design a suitable filter (e.g., use Remez approach of lab 4)
f0 = Fs/12;

% number of samples is N+1
N = 100; % must be an even number
% limit frequencies
al = 0.1; % transition bandwidth in percentage
fp = f0*(1-al); % pass band upper limit
fs = f0; % stop band lower limit
err_lim = 0.0037; % -48.636 dB
[NN,Fo,Ao,W] = firpmord([fp fs],[1 0],[1 1]*err_lim,Fs);
% define filter
h0 = firpm(NN,Fo,Ao,W);

% behaviour in frequency obtained by use of function freqz
[H0, ff]= freqz(h0,Nx,'whole',Fs);

% filter the signal
z = filter(h0,1,x);
% show the absolute value in frequency of the filtered signal
Z = T*fft(z);
figure(2);
plot(frequency/1e3,20*log10(abs(Z))); grid;
xlabel('frequency [kHz]'); xlim([0 Fs/1e3]);
title('original audio signal in frequency');


%% sample signal

% sample signal


% show absolute value in frequency of the sampled signal


% listen to first 10 seconds of the sampled signal
Fs1 = Fs/6; % this is its sampling rate

