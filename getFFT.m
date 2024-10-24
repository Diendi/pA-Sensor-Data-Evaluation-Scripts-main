%performs fft and returns result in input unit. can be manually converted
%to dB by user
function [f, data] = getFFT(Fs, U)

U=U(~isnan(U));

L=length(U);

Y1=fft(U);
P21 = abs((Y1/L));
P11 = P21(1:L/2+1);
P11(2:end-1) = 2*P11(2:end-1);

f = Fs*(0:(L/2))/L;
data = P11;

end