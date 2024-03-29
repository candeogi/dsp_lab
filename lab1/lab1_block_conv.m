clear all
close all
clc


%% Lab 1 BLOCK CONVOLUTION


%% 1. signal and filter

% load ECG data
load('data_ecg')
T = 1/125; % sampling period
x = ecg-mean(ecg); % ECG signal
tx = T*(0:length(x)-1); % time instants associated to x

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
y = exp(-0.25*ty/T); y = y/sum(y)/T; % exponential signal
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
Nx = length(x);
Ny = length(y);
N = Nx+Ny-1;
x1 = [x, zeros(1,N-Nx)]; % add zeros
y1 = [y, zeros(1,N-Ny)]; % add zeros
X1 = T*fft(x1);  % Fourier transform of x
Y1 = T*fft(y1);  % Fourier transform of y
Z1 = X1.*Y1;     % Fourier transform of z (.* is mult entry b entry)

z1 = ifft(Z1)/T; % convolution signal

hold on
plot(tz,z1,'r-.') % plot the signal
hold off
%legend('via conv','via fft')
%disp(['norm of the difference of the two approaches = ' num2str(norm(z-z1))])


%% 3. block convolution - overlap and add

%N = 4019
%Nx = 4000
%Ny = 20

%M = 20 blocchi lunghi 200
M = 10*Ny; % exactly matches with Nx = 20*M1
nOfBlocks = Nx/M;

%blocks is a 1x200 cell that contains all the blocks to be computed for the
%overlap and add method
blocks  = {};
%BLOCKS contains the blocks on freq domain
BLOCKS = {};

%for the FFT now N = M + Ny - 1
NN = M + Ny - 1;
padding = zeros(1,Ny-1);

newy1 =[y zeros(1, NN - Ny)];
newY1 = T*fft(newy1); %lungo 219

%i va da 1 a 20
for i = 1:nOfBlocks
    %0*20+1
    %1*20+1
    startIndexBlock = ((i-1)*M+1);
    endIndexBlock = (i*M);
    %pad the blocks in preparation for the fft
    blocks{1,i} = [x(startIndexBlock : endIndexBlock) padding]; %lungo 219
    %compute the fft on the padded blocks
    BLOCKS{1,i} = T*fft(blocks{1,i}); %theres no more zeros at the end
    %lets use the same structure for the product 
    BLOCKS{1,i} = BLOCKS{1,i}.*newY1;
    blocks{1,i} = ifft(BLOCKS{1,i})/T;
end
%now i need to add the overlapping parts
%each block is long NN = 219

result = zeros(1,Nx +Ny -1);

for i=1:nOfBlocks
    for j=1:length(blocks{1,i})
        % 0*200+1
        % 0*200+2 etc...
        % ...
        % 0*200+219 fine 
        result((i-1)*M+j) = result((i-1)*M+j) + blocks{1,i}(j);   
    end
end

hold on
plot(tz,result ,'g:') % plot the signal
hold off


%% 4. block convolution - overlap and save

%lets add the Ny-1 zeros ad the start and end for the transient
padded_x = [zeros(1,Ny-1) x zeros(1,Ny-1)];

%divide it in overlapping blocks of length N=M+Ny-1
blockLength = M+Ny-1;

blocks2 = {};
BLOCKS2= {};

padded_y = [y zeros(1,M-1)];
PADDED_Y = T*fft(padded_y);

result3=[];

for i=1:20
    startidx=(i-1)*M+1;
    endidx = startidx+(blockLength-1);
    %make the first block cut
    %blocks2{1,i}=padded_x(startidx : endidx);
    %alternative way to write this:
    blocks2{1,i}=padded_x((1:blockLength)+((i-1)*M));
    
    
    
    %now cyclically convolve with y then discard the first 19 samples
    %fft
    BLOCKS2{1,i}=T*fft(blocks{1,i});
    %multiplication with fft of y
    BLOCKS{1,i}=BLOCKS2{1,i}.*PADDED_Y;
    %back to time domain
    blocks2{1,i}=ifft(BLOCKS{1,i})/T;
    %discard the first Ny-1 values
    blocks2{1,i}=blocks2{1,i}(Ny:end);
    %end
    
    %result3 = [result3, blocks2{1,i}(Ny:end)];
end

%result2 = zeros(1,Nx +Ny -1);
result2 = [];
for i=1:20
    for j=1:length(blocks2{1,i}) 
        %result2((1-1)*200+1)
        %result2((i-1)*M+j) = blocks2{1,i}(j);   
        %result2 = [result2, blocks2{1,i}(j)];
    end
    result2 = [result2, blocks2{1,i}];
end

hold on
plot(tz,result2,'m--') % plot the signal
hold off
legend('via conv','via fft','via overlap','via o n s')
disp(['norm of the difference of the two approaches = ' num2str(norm(z-result2))])

%%  cylic convolution function via fft: x and y must have the same length
%io praticamente facevo sta roba a mano ogni volta

function z = cyclic_conv(x,y,T)
    X = T*fft(x);
    Y = T*fft(y);
    Z = X.*Y;
    z = ifft(Z)/T;
end

