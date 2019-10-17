
%% bluestein's fft function %%%%%%%%%%%%%%%

%s are the samples
%f0 is the starting point
%dF sampling period in the frequency domain
%T i guess is the sampl period. 1/125
%
%returns Sk:samples of the fft, fk: frequencies.
function [Sk, fk] = bluestein(s,f0,dF,T)

    N = length(s); % number of samples
    fk = f0+dF*(0:N-1); % frequency samples of interest
    n = T*(0:length(s)-1); % time instants associated to s
    a=dF*T; %parameter for chirp signal for the bluestein approach
    
    %multiply by exp function

    z = s.*exp(-1i.*(((2*pi*f0*T).*n) + ((pi*a).*n.^2)));
    
    z = [z zeros(1,N)]; %add N zeros - this is vector q
    Z = T*fft(z);
    
    n = T*(-N:N-1);
    %generate vector q
    q = exp((1i*pi*a).*(n).^2);
    Q = T*fft(q);
    
    %correct by multiplying with this exp
    Sk = Q.*Z;
    Sk = Sk(N+1:end); %discard the first N samples

    
end