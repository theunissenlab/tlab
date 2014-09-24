function f = fftfreqs(n,fs)

%% Copy of numpy.fft.helper.fftfreq()
% n is the number of points
% fs is the sample frequency

switch mod(n,2)
    case 0% Even window length
        n1 = (n/2) - 1;
        n2 = -n/2;
    case 1% Odd window length
        n1 = (n-1)/2;
        n2 = (1-n)/2;
end

f = [0:n1,n2:-1] * fs/n;