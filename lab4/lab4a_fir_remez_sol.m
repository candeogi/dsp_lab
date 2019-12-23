clearvars
close all
clc

%% LAB 4 FIR FILTER - LOW PASS - REMEZ ALGORITHM TECHNIQUE
%
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

% firpmord  Parks-McClellan optimal equiripple FIR order estimator.
% 
%     [N,Fo,Ao,W] = firpmord(F,A,DEV,Fs) 
%
%     finds the approximate order N, 
%     normalized frequency band edges Fo, frequency band magnitudes Ao and 
%     weights W to be used by the FIRPM function as follows:
%
%         h = firpm(N,Fo,Ao,W)
% 
%     The resulting filter will approximately meet the specifications given
%     by the input parameters F, A, and DEV.  F is a vector of cutoff
%     frequencies in Hz, in ascending order between 0 and half the sampling
%     frequency Fs. If you do not specify Fs, it defaults to 2.  A is a
%     vector specifying the desired function's amplitude on the bands defined
%     by F. The length of F is twice the length of A, minus 2 (it must
%     therefore be even).  The first frequency band always starts at zero,
%     and the last always ends at Fs/2.  It is not necessary to add these
%     elements to the  F vector.  DEV is a vector of maximum deviations or
%     ripples (in linear units) allowable for each band.  DEV must have the
%     same length as A.

err_lim = 0.0037; % -48.636 dB
[NN,Fo,Ao,W] = firpmord([fp fs],[1 0],[1 1]*err_lim,Fp);
disp(['firpmord suggests a filter of order ' num2str(NN) ...
      ' for guaranteeing a ' num2str(20*log10(err_lim)) ' dB error'])

% define filter
h0 = firpm(N,Fo,Ao,W)/T;
t = T*(-N/2:N/2);

% behaviour in frequency obtained by use of function freqz
[H0,ff] = freqz(h0,1,8*(N+1),Fp);
H0 = T*H0; % normalization factor

% show results
figure(1)
subplot(2,1,1)
stem(t,h0); grid; title('FIR - remez approach - time domain')
subplot(2,1,2)
plot(ff,20*log10(abs(H0))); grid; xlim([0 Fp/2]); ylim([-60, 5]); 
title('frequency domain')
hold on; plot([1,1]*fp,ylim,'r--'); plot([1,1]*fs,ylim,'r--'); hold off;

