clearvars
close all
clc

%% LAB 5 IIR FILTER - LOW PASS - CHEBYSHEV 1 FILTER TECHNIQUE

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
Rp = -20*log10(1-2/A); % [dB] ripple in the pass band (from FIR)
ep = sqrt(10^(Rp/10)-1); % ripple constant epsilon
disp(['Rp = ' num2str(Rp) ' dB, Rs = ' num2str(Rs) ' dB'])


%%  build Chebyshev filter coefficients 

% define constants
fp_1 = tan(pi*fp*T)/2*pi;
fs_1 = tan(pi*fs*T)/2*pi;
N = log((ep^2)/(A^2-1))/(2*log(fp_1/fs_1));
omega0 = 2*pi*fp_1/(ep^(1/N));

% estimate poles in the s domain
sk = [];
for k = 0:(N-1)
   s 
end

% map poles in the z domain


% graphically show the location of poles and zeros


% show the resulting frequency response
%
% suggestion: evaluate the frequency responses of single layers (made by
% 1 zero + 1 pole) over the frequency range, and then multiply them to
% observe the overall behavior of the filter over frequencies




%%  compare with MatLab built-in solution 
%   you can also use function: filterDesigner

% % estimate the order
% [NCI, wp] = cheb1ord(fp/(Fp/2), fs/(Fp/2), Rp, Rs);
% disp(['The order estimated by MatLab is N = ' num2str(NCI)])
% % extract coefficients
% [B, A] = cheby1(NCI, Rp, wp);
% fvtool(B, A);

