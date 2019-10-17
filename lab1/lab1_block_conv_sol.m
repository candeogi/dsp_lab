clear all
close all
clc


%% Lab 1 BLOCK CONVOLUTION


%% 1. signal and filter

% load ECG data
load('data_ecg')
T = 1/125; % sampling period
testmean = mean(ecg)
x = ecg-mean(ecg); % ECG signal
teststuff = (0:length(x)-1);
tx = T*(0:length(x)-1); % time instants associated to x
Nx = length(x);

% plot the ECG signal in the time domain
figure(1)
subplot(2,2,1)
plot(tx,x) % plot the signal
grid % activate the grid
xlim([0,20]) % zoom on a signal portion
xlabel('time [s]'); ylabel('signal x');
title('The ECG signal in the time domain')

% filter
ty = (0:19)*T; 
y = exp(-0.25*ty/T); 
y = y/sum(y)/T; % exponential signal
Ny = length(y);
subplot(2,2,2)
plot(ty,y) % plot the signal
grid % activate the grid
xlabel('time [s]'); ylabel('signal y');
title('The low-pass filter')


%% 2. convolution 

% in the time domain
z = T*conv(x,y); % discrete convolution
tz = tx(1)+ty(1):T:tx(end)+ty(end); % property of convolution
subplot(2,1,2)
plot(tz,z) % plot the signal
grid % activate the grid
xlim([0,20]) % zoom on a signal portion
xlabel('time [s]'); ylabel('signal z');
title('The convolution')

% in the frequency domain
x1 = [x, zeros(1,Ny-1)]; % add zeros
y1 = [y, zeros(1,Nx-1)]; % add zeros
z1 = cyclic_conv(x1,y1,T); % cyclic convolution via fft 
                           % (see the function definition below)
hold on
plot(tz,z1,'r-.') % plot the signal
hold off
legend('via conv','via fft')
disp(['difference of conv versus fft = ' num2str(norm(z-z1))])


%% 3. block convolution - overlap and add

M = 10*Ny; % block length; it provides 20 blocks since Nx = 20*M

N = M+Ny-1; % block constant
z2 = zeros(1,Nx+Ny-1); % prepare an all-zero vector for the output

for i = 1:Nx/M % cycle per each block
    b1 = [x((1:M)+(i-1)*M), zeros(1,Ny-1)]; % build block and add zeros
    y1 = [y, zeros(1,M-1)]; % add zeros to y
    c1 = cyclic_conv(b1,y1,T); % cyclic convolution via fft 
                               % (see the function definition below)
    pos = (i-1)*M+(1:N); % positions of the output vector to update
    z2(pos) = z2(pos) + c1; % add in the correct place
end

disp(['difference of conv versus overlap-and-add = ' num2str(norm(z-z2))])


%% 4. block convolution - overlap and save

x2 = [zeros(1,Ny-1), x, zeros(1,N)]; % add zeros before and after
z3 = []; % prepare an empty vector for the output

for i = 1:length(x2)/M % cycle per each block
    b1 = x2((1:N)+(i-1)*M); % build block, do not add zeros
    y1 = [y, zeros(1,M-1)]; % add zeros to y
    c1 = cyclic_conv(b1,y1,T); % cyclic convolution via fft 
                               % (see the function definition below)
    z3 = [z3, c1(Ny:end)]; % append the new output, keep only part of the samples
end

z3 = z3(1:Nx+Ny-1); % restrict the vector to the active samples
disp(['difference of conv versus overlap-and-save = ' num2str(norm(z-z3))])


%%  cylic convolution function via fft: x and y must have the same length

function z = cyclic_conv(x,y,T)
    X = T*fft(x);
    Y = T*fft(y);
    Z = X.*Y;
    z = ifft(Z)/T;
end


