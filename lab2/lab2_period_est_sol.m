clearvars
close all
clc


%% Lab 2 - ESTIMATE THE ECG PERIOD


%% 1. time domain 

% load ECG data
load('data_ecg')
s = ecg-mean(ecg); % ECG signal
T = 1/125; % sampling period
N = length(s); % length of the signal
t = T*(0:N-1); % time instants associated to s

% empirical detection method in time: find the distance between two maxima, 
% divide by the number of periods between them
[m1,pos] = max(s(t<2)); t1 = t(t<2); t1 = t1(pos);
[m2,pos] = max(s(t>30)); t2 = t(t>30); t2 = t2(pos);
disp(['period estimate in time Tp = ' num2str((t2-t1)/23)])

% plot the ECG signal in the time domain
figure
subplot(2,1,1)
plot(t,s) % plot the signal
hold on; plot([t1 t2],[m1 m2],'rx'); hold off; 
grid % activate the grid
xlim([0,max(t)]) % zoom on a signal portion
title(['Period estimate in the time domain: Tp = ' num2str((t2-t1)/23)])


%% 2. frequency domain 

% apply windowing
window = 'Kaiser';
switch window
    case 'Hann'
        w = 0.5+0.5*cos(2*pi*(-N/2:N/2-1)/N);
    case 'Hamming'
        w = 0.54+0.46*cos(2*pi*(-N/2:N/2-1)/N);
    case 'Blackman'
        w = 0.42+0.5*cos(2*pi*(-N/2:N/2-1)/N)+0.08*cos(4*pi*(-N/2:N/2-1)/N);
    case 'Kaiser'
        be = .1102*(50-8.7);
        w = besseli(0,be*sqrt(1-(2*(-N/2:N/2-1)/N).^2))/besseli(0,be); % Kaiser
end
sw = [s.*w, zeros(1,2^ceil(log2(N))-N)]; % windowed + zero padding 
s = [s, zeros(1,2^ceil(log2(N))-N)]; % zero padding up to a power of 2

% extract frequency samples by Bluestein's algorithm
dF = 1/(T*N*256); % frequency samples 
f0 = 0.5; % starting point
[Sk, fk] = bluestein(s,f0,dF,T);
[Swk, fk] = bluestein(sw,f0,dF,T);

% empirical detection method in frequency: find the frequency of the
% maximum, its inverse is the period
[S1,pos] = max(abs(Sk)); f1 = fk(pos);
[S2,pos] = max(abs(Swk)); f2 = fk(pos);
disp(['period estimate in frequency (original) Tp = ' num2str(1/f1)])
disp(['period estimate in frequency (windowed) Tp = ' num2str(1/f2)])

% plot the ECG signal in the frequency domain
subplot(2,2,3)
semilogy(fk,abs(Sk)) % plot the Fourier transform
hold on; semilogy(f1,S1,'rx'); hold off;
grid % activate the grid
title(['in frequency (original): Tp = ' num2str(1/f1)])
subplot(2,2,4)
semilogy(fk,abs(Swk)) % plot the Fourier transform
hold on; semilogy(f2,S2,'rx'); hold off;
grid % activate the grid
title(['in frequency (' window ' window): Tp = ' num2str(1/f2)])


%% functions %%%%%%%%%%%%%%%

function [Sk, fk] = bluestein(s,f0,dF,T)

N = length(s); % number of samples
fk = f0+dF*(0:N-1); % frequency samples of interest
% intermediate signal z
z = s.*exp(-2i*pi*f0*T*(0:N-1)-1i*pi*T*dF*(0:N-1).^2); 
z = [z, zeros(1,N)];
% intermediate signal q
q = exp(1i*pi*T*dF*(-N:N-1).^2); 
% cyclic convolution via FFT
v = T*ifft(fft(z).*fft(q));
% extract a part and correct
Sk = v(N+1:end).*exp(-1i*pi*T*dF*(0:N-1).^2);
end