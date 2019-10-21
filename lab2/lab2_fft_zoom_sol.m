clearvars
close all
clc


%% Lab 2 FFT ZOOM


%% 1. standard FFT

% load ECG data
load('data_ecg')
s = ecg-mean(ecg); % ECG signal
T = 1/125; % sampling period

% zero padding up to a power of 2
N = length(s); % length of the signal
s = [s, zeros(1,2^ceil(log2(N))-N)]; % zero padding up to a power of 2
t = T*(0:length(s)-1); % time instants associated to s

% plot the ECG signal in the time domain
figure(1)
plot(t,s) % plot the signal
grid % activate the grid
xlim([0,20]) % zoom on a signal portion
title('The ECG signal in the time domain')

% extract frequency samples
S = T*fft(s); % apply the DFT operator
N = length(S); % number of samples
F = 1/(T*N); % sampling period in the frequency domain
f = F*(0:N-1); % frequency instants associated to S

% plot the ECG signal in the frequency domain
figure(2)
semilogy(f,abs(S)) % plot the Fourier transform
grid % activate the grid
ylim([1e-2,1e4]) % zoom vertically
title('The ECG signal in the frequency domain')

% zoom on lower frequencies
figure(3)
semilogy(f,abs(S),'.-') % plot the Fourier transform
grid % activate the grid
ylim([1e1,2e3]) % zoom vertically
xlim([0, 5]) % zoom horizontally
title('The ECG signal in the frequency domain (zoom)')


%% 2. Interpolation in frequency by zero padding

% extract frequency samples
disp('Running zero-padding fft ...')
tic % starts the counter
s1 = [s zeros(1,63*length(s))]; % zero padded signal
S = T*fft(s1); % apply the DFT operator
toc % reads the counter
N = length(S); % number of samples
F = 1/(T*N); % sampling period in the frequency domain
f = F*(0:N-1); % frequency instants associated to S

% plot the result along with the standard one
figure(3)
hold on
semilogy(f,abs(S)) % plot the Fourier transform
hold off


%% 3. Zooming by Bluestein's algorithm

% determine the frequency samples of interest
dF = F; % same frequency samples as per the zero padded solution
f0 = 0*F; % starting point

% run Bluestein's algorithm
disp('Running Bluestein''s algorithm...')
tic % starts the counter
[Sk, fk] = bluestein(s,f0,dF,T);
toc % reads the counter

% plot the result along with the standard one
figure(3)
hold on
semilogy(fk,abs(Sk)) % plot the Fourier transform
hold off
legend('original signal','zero padded','Bluestein''s algorithm')

disp(['norm of the difference between Bluestein and zero-padding is ' num2str(norm(Sk-S(1:length(Sk))))])



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