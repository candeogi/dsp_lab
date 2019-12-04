clearvars
close all
clc

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
F = (Fp/(N+1))/32; % min 32 samples per cosine period - 32 è empirico, 
% N+1 sarebbe sbagliato, sarebbe N e basta
f = [0:F:fp, fs:F:Fp/2].'; % frequency samples, column vector
f_len = length(f);

% 1) build matrices for the linear program function

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

% 2) use MatLab function x = linprog(g,A,b) to solve the 
% linear programming problem:
% 
%             min g'*x    subject to:   A*x <= b
%              x
%
% where x(1:N/2+1) carries the filter samples, and
% x(N/2+2) is the error value delta

result_pl = linprog(obj_f,A,b);
half_h0 = result_pl(1:(end-1)); %51 samples , stem(x) to see
delta = result_pl(end);

% 3) show results in time (h0) and frequency (H0)

figure(1)

flipped_half_h0 = flipud(half_h0);
h0 = [half_h0;flipped_half_h0(2:end)];
stem(h0);
t = [-(N/2)*T:T:(N/2)*T];

%top plot
subplot(2,1,1);
stem(t,h0); grid;
xlabel('time [t]'); ylabel('filter h0');
title('filter in time domain');

%X = V*half_h0;
%from solutions this func is used
[H0,ff] = freqz(h0,1,8*(N+1),Fp);

%bottom plot
subplot(2,1,2);
plot(ff,20*log10(abs(H0)));
grid;
xlabel('frequency [f]'); ylabel('filter H0');
title('filter in frequency domain')

