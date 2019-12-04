clearvars
close all
clc

%% LAB 4 FIR FILTER - LOW PASS - BUTTERWORTH FILTER TECHNIQUE
%

% ideally a low-pass filter with cutoff frequency 
f0 = 1e3; % [Hz]
% and repetition rate
Fp = 8e3; % [Hz]
% so that the samplig period is
T = 1/Fp; % [s]

% limit frequencies
al = 0.1; % transition bandwidth in percentage
fp = f0*(1-al); % pass band upper limit
fs = f0*(1+al); % stop band lower limit

% limits (the same as per the FIR approach)
Rs = 48.636; % [dB] attenuation in the stop band (from FIR)
A = 10^(Rs/20); % attenuation constant A
Rp = -20*log10(1-1/A); % [dB] ripple in the pass band (from FIR)
ep = sqrt(10^(Rp/10)-1); % ripple constant epsilon
disp(['Rp = ' num2str(Rp) ' dB, Rs = ' num2str(Rs) ' dB'])

% map back into s-domain constraints
fp1 = tan(pi*fp*T)/(2*pi);
fs1 = tan(pi*fs*T)/(2*pi);

% estimate the order N
N = ceil(0.5*log10((A^2-1)/ep^2)/log10(fs1/fp1)); 
disp(['The estimated order is N = ' num2str(N)])

% estimate the 3dB frequency
f01 = fp1/(ep^(1/N));

% poles
sk = 2*pi*f01*exp(1i*pi/2*(1+(1:2:2*N-1)/N));
pk = (1+sk)./(1-sk);

% builds the frequency response
f = 0:Fp/1000:Fp;  % 
H = ones(size(f));
for k = 1:N
    num = (1+exp(-2i*pi*f*T)) / 2; % normalized at f=0
    den = (1-pk(k)*exp(-2i*pi*f*T)) / (1-pk(k)); % normalized at f=0
    H = H.*num./den;
end

% show poles and zeros
figure(1)
zplane([1 1],[1]); % Matlab notation for the denominator 
hold on; plot(pk,'x'); hold off

% show frequency plot
figure(2); 
plot(f,20*log10(abs(H)));
grid on; axis([0 Fp -60 5]); xlabel('f [Hz]'); ylabel('|H| [dB]'); 
hold on; plot([1,1]*fp,ylim,'r--'); plot([1,1]*fs,ylim,'r--'); 
plot(xlim,-[1,1]*Rs,'r--'); plot(xlim,-[1,1]*Rp,'r--'); hold off;


%%  compare with MatLab built-in solution 
%   you can also use function: filterDesigner

% estimate the order
[NB, w0] = buttord(fp/(Fp/2), fs/(Fp/2), Rp, Rs);
disp(['The order estimated by MatLab is N = ' num2str(N)])
% extract coefficients
[B, A] = butter(NB, w0);
fvtool(B, A);
% plot result
[HB, W] = freqz(B, A, 1024*2);
figure(2); hold on; plot(W/(2*pi*T),20*log10(abs(HB))); hold off
