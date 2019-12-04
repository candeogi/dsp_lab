clearvars
close all
clc

%% LAB 4 FIR FILTER - BAND PASS


% Repetition rate
Fp = 8e3; % [Hz]
% so that the samplig period is
T = 1/Fp; % [s]

% pass band from fp1 to fp2
fp1 = 500; % [Hz] from
fp2 = 1500; % [Hz] to

% attenuation bands
fs1 = 300; % [Hz] stopband #1: from 0 to fs1
fs2 = 1700; % [Hz] stopband #2: from fs2 to Fp/2


% 1) use REMEZ ALGORITHM: firpmor and firpm
err_lim = [1e-2 1e-3 1e-2]; % required limits (delta)

F=[fs1,fp1,fp2,fs2];
A=[0,1,0];
DEV=err_lim;
[N,Fo,Ao,W] = firpmord(F,A,DEV,Fp);
h = firpm(N,Fo,Ao,W)/T;


% illustrate filter in time and frequency
figure
subplot(2,1,1)
t = T*(-N/2:N/2);
stem(t,h); grid; 
ylim([-2000 3000]);
title('FIR - linear program approach - time domain')
% behaviour in frequency obtained by use of function freqz
[H0,ff] = freqz(h,1,8*(N+1),Fp);
H0 = T*H0; % normalization factor
subplot(2,1,2)
plot(ff,20*log10(abs(H0))); grid; xlim([0 Fp/2]);
ylim([-80 20]); 
title('frequency domain')

% 2) use the LINEAR PROGRAM approach: linprog

%reset N
N = 100;
F = Fp/(N/2)/32; % min 32 samles per cosine period
f = [0:F:fs1, fp1:F:fp2, fs2:F:Fp/2].'; % frequency samples, column vector
f_len = length(f);

%matrix V
V = [];
for k=N/2:-1:1
    cosColumn = 2*T*cos(2*pi*f*k*T);
    V = [V cosColumn];
end
V = [V T*ones(f_len,1)];

%vector r
first = [0:F:fs1];
lenFirst=length(first);
second = [fp1:F:fp2];
lenSecond=length(second);
third= [fs2:F:Fp/2];
lenThird= length(third);
r_first = zeros(lenFirst,1);
r_second = ones(lenSecond,1);
r_third = zeros(lenThird,1);
r = [r_first;r_second;r_third];
r = transpose(r);

%vector w
w=ones(f_len,1);

%build matrixes for linprog function
A = [-V -w; V -w];
b = [-r r];

%g in the linprog function
g = [zeros(ceil(N/2)+1,1);1];
x = linprog(g,A,b);

% illustrate filter in time and frequency
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