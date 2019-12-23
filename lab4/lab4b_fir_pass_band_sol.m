clearvars
close all
clc

%% LAB 4 FIR FILTER - BAND PASS


% Repetition rate
Fp = 8e3; % [Hz]
% so that the samplig period is
T = 1/Fp; % [s]

% pass band
fp1 = 500; % [Hz] from
fp2 = 1500; % [Hz] to

% attenuation bands
fs1 = 300; % [Hz] from 0 to
fs2 = 1700; % [Hz] from ... to Fp/2

% required limits
err_lim = [1e-3 1e-2 1e-3];


%% 1) REMEZ ALGORITHM TECHNIQUE

% define filter
[N,Fo,Ao,W] = firpmord([fs1 fp1 fp2 fs2],[0 1 0],err_lim,Fp); % identify order
N = 2*ceil(N/2); % must be an even number; number of samples is N+1
disp(['firpmord suggests a filter of order ' num2str(N)])
h0 = firpm(N,Fo,Ao,W)/T; 
t = T*(-N/2:N/2);

% behaviour in frequency obtained by use of function freqz
[H0,ff] = freqz(h0,1,8*(N+1),Fp);
H0 = T*H0; % normalization factor

% show results
figure(1)
subplot(2,1,1)
stem(t,h0); grid; title('FIR - high pass - Remez approach - time domain')
subplot(2,1,2)
plot(ff,20*log10(abs(H0))); grid; xlim([0 Fp/2]); ylim([-80, 5]); 
title('frequency domain')
hold on; plot([1,1]*fp1,ylim,'r--'); plot([1,1]*fp2,ylim,'r--'); 
plot([1,1]*fs1,ylim,'r--'); plot([1,1]*fs2,ylim,'r--'); hold off;


%% 2) LINEAR PROGRAM TECHNIQUE

% frequencies samples of interest
F = Fp/(N+1)/32; % min 32 samles per cosine period
f = [0:F:fs1, fp1:fp2, fs2:F:Fp/2].'; % frequency samples, column vector

% build matrices
d = (f>=fp1).*(f<=fp2); % ideal filter shape
w = err_lim(1)*(f<=fs1) + err_lim(2)*(f>=fp1).*(f<=fp2) ...
    + err_lim(3)*(f>=fs2); % weighting function
V = T*ones(size(f)); % cosines matrix
for n = 1:N/2
    V = [V,2*T*cos(2*pi*f*n*T)];
end

% linear programming solution
g = [zeros(N/2+1,1);1];
A = [-V, -w; V, -w];
b = [-d;d];
x = linprog(g,A,b); % solve the linear program

% define filter
h0 = [x(N/2+1:-1:2);x(1:N/2+1)];
t = T*(-N/2:N/2);

% behaviour in frequency obtained by use of function freqz
[H0,ff] = freqz(h0,1,8*(N+1),Fp);
H0 = T*H0; % normalization factor

% show results
figure(2)
subplot(2,1,1)
stem(t,h0); grid; title('FIR - high pass - linear programming approach - time domain')
subplot(2,1,2)
plot(ff,20*log10(abs(H0))); grid; xlim([0 Fp/2]); ylim([-80, 5]); 
title('frequency domain')
hold on; plot([1,1]*fp1,ylim,'r--'); plot([1,1]*fp2,ylim,'r--'); 
plot([1,1]*fs1,ylim,'r--'); plot([1,1]*fs2,ylim,'r--'); hold off;
