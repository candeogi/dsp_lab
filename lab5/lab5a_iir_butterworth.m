clearvars
close all
clc

%% LAB 5 IIR FILTER - LOW PASS - BUTTERWORTH FILTER TECHNIQUE
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


%%  build Butterworth filter coefficients 

% define constants
fp_1 = tan(pi*fp*T)/(2*pi);
fs_1 = tan(pi*fs*T)/(2*pi);

%N = ceil(log10((ep^2)/(A^(2)-1))/(2*log10(fp_1/fs_1)));
N = ceil(0.5*log10((A^2-1)/ep^2)/log10(fs_1/fp_1)); 
disp(['The estimated order is N = ' num2str(N)])
omega0 = (2*pi*fp_1)/(ep^(1/N));

% estimate poles in the s domain
sk = 1i*omega0*exp(1i*pi*(1+2*(0:1:(N-1)))/(2*N));

% map poles in the z domain
pk = [];
for k = 1:N
    pk = [pk , (1+sk(k))/(1-sk(k))];
end 

% graphically show the location of poles and zeros
figure(1);
zplane([1 1],[1]);
hold on;
plot(pk,'x'); 
hold off;

% show the resulting frequency response
f = 0:Fp/1000:Fp; %scelti a caso i guess
H = ones(size(f));
for k = 1:N
    num = (1+exp(-2i*pi*f*T)) / 2; % normalized at f=0
    den = (1-pk(k)*exp(-2i*pi*f*T)) / (1-pk(k)); % normalized at f=0
    H = H.*num./den;
end

% suggestion: evaluate the frequency responses of single layers (made by
% 1 zero + 1 pole) over the frequency range, and then multiply them to
% observe the overall behavior of the filter over frequencies

figure(2);
plot(f,20*log10(abs(H)));
grid on;
axis([0 Fp -60 5]);
xlabel('f [Hz]');
ylabel('|H| [dB]');
hold on; 
plot([1,1]*fp,ylim,'r--');
plot([1,1]*fs,ylim,'r--'); 
plot(xlim,-[1,1]*Rs,'r--'); 
plot(xlim,-[1,1]*Rp,'r--'); 
hold off;



%% compare with MatLab built-in solution (uncomment the following)
%   you can also use function: filterDesigner

% % estimate the order
% [NB, w0] = buttord(fp/(Fp/2), fs/(Fp/2), Rp, Rs);
% disp(['The order estimated by MatLab is N = ' num2str(N)])
% % extract coefficients
% [B, A] = butter(NB, w0);
% % graphycally show the result
% fvtool(B, A);
