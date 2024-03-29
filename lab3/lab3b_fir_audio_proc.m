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

% show the ORIGINAL signal in time and frequency
Ny = length(y); % length
ty = T*(0:Ny-1); % time samples
Y = T*fft(y); % fft
fy = (0:Ny-1)/(T*Ny); % frequency samples
figure(1)
subplot(2,1,1) % show a portion of time
plot(ty,y); grid; xlim([0.1 0.2]); xlabel('time [s]'); 
title('original audio signal in time')
subplot(2,1,2) % show frequency content in dB scale
plot(fy/1e3,20*log10(abs(Y))); grid; xlim([0 Fp/2e3]); ylim([-150 -20]) 
xlabel('frequency [kHz]'); title('original audio signal in frequency')

% listen to the original content
%disp('press any key to listen to the ORIGINAL audio signal')
%pause
%sound(y(1:5*Fp), Fp); % play the first 5 seconds 
%pause(5)

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
F = Fp/(N+1)/32; % min 32 samples per cosine period
f = [0:F:fp, fs:F:Fp/2].'; % frequency samples, column vector

% 1) build matrices for the linear program function
f_len = length(f);


%matrix V
V = [];
for k=N/2:-1:1
    cosColumn = 2*T*cos(2*pi*f*k*T);
    V = [V cosColumn];
end
V = [V T*ones(f_len,1)];

%vector r
first = [0:F:fp];
lenFirst=length(first);
second = [fs:F:Fp/2];
lenSecond=length(second);
r_first = ones(lenFirst,1);
r_second = zeros(lenSecond,1);
r = [r_first;r_second];
r = transpose(r);
%vector w
w = ones(f_len,1);
%w=transpose(w);

%build matrixes for linprog function
A = [-V -w; V -w];
b = [-1*r r];

%g in the linprog function
obj_f = [zeros(1,N/2+1) 1];

% 2) use linear programming technique
h_delta = linprog(obj_f,A,b);
h0 = h_delta(1:(end-1)); %51 samples , stem(x) to see
delta = h_delta(end);

% 3) show results 
figure(2)
flipped_h = flipud(h0);
h_mirror = [h0;flipped_h(2:end)];
stem(h_mirror);
time_axis = [-(N/2)*T:T:(N/2)*T];

%top plot
subplot(2,1,1);
stem(time_axis,h_mirror);
grid;
xlabel('time [t]'); ylabel('filter h0');
title('filter in time domain');

%bottom plot
subplot(2,1,2);
H = V*h0;
plot(f,H);
grid;
xlabel('frequency [f]'); ylabel('filter H0');
title('filter in frequency domain')

%% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%% use the filter on the audio signal

% 1) filter the signal
z = T*conv(y,h0);

% show the FILTERED signal in time and frequency
Nz = length(z);
tz = T*(0:Nz-1); % time samples
Z = T*fft(z); % fft
fz = (0:Nz-1)/(T*Nz); % frequency samples
figure(3)
subplot(2,1,1)
plot(tz,z); grid; xlim([0.1 0.2]); xlabel('time [s]'); 
title('filtered audio signal in time')
subplot(2,1,2)
plot(fz/1e3,20*log10(abs(Z))); grid; xlim([0 Fp/2e3]); ylim([-150 -20])
xlabel('frequency [kHz]'); title('filtered audio signal in frequency')
% 2) show the FILTERED signal in time and frequency


% 3) listen to the filtered content 
