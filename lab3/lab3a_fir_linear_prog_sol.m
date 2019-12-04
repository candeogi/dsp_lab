clearvars
close all
clc

%% LAB 3 FIR FILTER - LOW PASS - LINEAR PROGRAMMING TECHNIQUE

%% PART 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% filter constants

% ideally a low-pass filter with cutoff frequency 
f0 = 1e3; % [Hz]
% and repetition rate
Fp = 8e3; % [Hz] 
% so that the samplig period is
T = 1/Fp; % [s]

% number of samples is N+1
N = 100; % must be an even number

% limit frequencies
al = 0.1; % transition bandwidth in percentage
fp = f0*(1-al); % pass band upper limit
fs = f0*(1+al); % stop band lower limit


%% PART 1 %%%%%%%%%%%%%%%%%%%%%%% linear programming design technique

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
        
% x = linprog(g,A,b) attempts to solve the linear
% programming problem:
% 
%             min g'*x    subject to:   A*x <= b
%              x
%
% where x(1:N/2+1) carries the filter samples, and
% x(N/2+2) is the error value delta
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
