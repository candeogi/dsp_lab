
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
    n =(0:N-1); % time instants associated to s
    a=dF*T; %parameter for chirp signal for the bluestein approach

    %multiply by exp function
    z = s.*exp(-2i*pi*f0*T*n-1i*pi*a*n.^2);
    z = [z, zeros(1,N)]; %add N zeros - this is vector q

    %generate vector q
    n2 =(-N:N-1);
    q = exp(1i*pi*a*n2.^2);

    %cyclically convolve them via fft
    v = T*ifft(fft(z).*fft(q));
    %discard the first N samples
    Sk = v(N+1:end); 
    %correct by multiplying with exp((-1i*pi*a).*fk.^2 )
    Sk = Sk.*(exp((-1i*pi*a)*n.^2 ));

end