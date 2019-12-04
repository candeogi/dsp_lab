clearvars
close all
clc

%% LAB 3 FIR FILTER - LOW PASS - LINEAR PROGRAMMING TECHNIQUE

%% PART 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% wave file & constants

% read a WAVE file (*.wav)
% y array containing sound samples y(n), n=1:length(y)
% Fp sampling frequency in Hz
[y, Fp] = audioread('chopin_pollini2.wav');
% so that the samplig period is
T = 1/Fp; % [s]

% show ORIGINAL signal in time and frequency
Ny = length(y); % length
ty = T*(0:Ny-1); % time samples
Y = T*fft(y); % fft
fy = (0:Ny-1)/(T*Ny); % frequency samples
figure
subplot(2,1,1) % show a portion of time
plot(ty,y); grid; xlim([0.1 0.2]); xlabel('time [s]'); 
title('original audio signal in time')
subplot(2,1,2) % show frequency content in dB scale
plot(fy/1e3,20*log10(abs(Y))); grid; xlim([0 Fp/2e3]); ylim([-150 -20]) 
xlabel('frequency [kHz]'); title('original audio signal in frequency')

% listen to the original content
disp('press any key to listen to the ORIGINAL audio signal')
pause
sound(y(1:5*Fp), Fp); % play the first 5 seconds 
pause(5)

% ideally apply a low-pass filter with cutoff frequency 
f0 = 2e3; % [Hz]
% limit frequencies
al = 0.1; % transition bandwidth in percentage
fp = f0*(1-al); % pass band upper limit
fs = f0*(1+al); % stop band lower limit
% number of samples is N+1
N = 100; % must be an even number


%% PART 1 %%%%%%%%%%%%%%%%%%%%%%% design your filter by linear programming

% frequencies samples of interest
F = Fp/(N/2)/32; % min 32 samples per cosine period
f = [0:F:fp, fs:F:Fp/2].'; % frequency samples, column vector

% build matrices for the linear program function
rect = @(x) 1*(abs(x)<0.5) + 0.5*(abs(x)==0.5); % inline function definition
d = rect(f/(2*f0)); % ideal filter shape
V = T*ones(size(f)); % cosines matrix
for n = 1:N/2
    V = [V,2*T*cos(2*pi*f*n*T)];
end
        
% use linear programming technique
g = [zeros(N/2+1,1);1];
A = [-V, -ones(size(f)); V, -ones(size(f))];
b = [-d;d];
x = linprog(g,A,b); % solve the linear program

% define filter
h0 = [x(N/2+1:-1:2);x(1:N/2+1)];
t = T*(-N/2:N/2);

% behaviour in frequency obtained by use of function freqz
[H0,ff] = freqz(h0,1,8*(N+1),Fp);
H0 = T*H0; % normalization factor

% show results
figure
subplot(2,1,1)
stem(t,h0); grid; title('FIR - linear program approach - time domain')
subplot(2,1,2)
plot(ff,20*log10(abs(H0))); grid; xlim([0 Fp/2]); ylim([-80, 5]); 
title('frequency domain')
hold on; plot([1,1]*fp,ylim,'r--'); plot([1,1]*fs,ylim,'r--'); 
plot(xlim,20*log10(x(end))*[1,1],'r--'); hold off;


%% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%% use the filter on the audio signal

% filter signal
z = T*conv(y,h0);

% show the FILTERED signal in time and frequency
Nz = length(z);
tz = T*(0:Nz-1); % time samples
Z = T*fft(z); % fft
fz = (0:Nz-1)/(T*Nz); % frequency samples
figure
subplot(2,1,1)
plot(tz,z); grid; xlim([0.1 0.2]); xlabel('time [s]'); 
title('filtered audio signal in time')
subplot(2,1,2)
plot(fz/1e3,20*log10(abs(Z))); grid; xlim([0 Fp/2e3]); ylim([-150 -20])
xlabel('frequency [kHz]'); title('filtered audio signal in frequency')

% listen to the filtered content
disp('press any key to listen to the FILTERED audio signal')
pause
sound(z(1:5*Fp), Fp); % play the first 5 seconds 
