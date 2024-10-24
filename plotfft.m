% plotfft.m
% Does the fft and plots it
% input is sampling frequency Fs and signal U, can also provide plot
% formatting
% Returns FFT result (f, FFT(U))
% Matthias Musil

%optionally enable windows here

function fftdata = plotfft(Fs,U,linlog,varargin);

%Fs=1/(t(2)-t(1));

%t=t(~isnan(U));
U=U(~isnan(U));

L=length(U);

Y1=fft(U);
P21 = abs(Y1/L);
P11 = P21(1:L/2+1);
P11(2:end-1) = 2*P11(2:end-1);
P11=20*log10(P11);

f = Fs*(0:(L/2))/L;
fftdata(:,1)=f;
fftdata(:,2)=P11;

cond = exist('linlog');

if cond
    if linlog == 'lin'
        plot(f,P11,varargin{:});
    end 
    if linlog == 'log'
        semilogx(f,P11,varargin{:});
    end
else
    plot(f,P11,varargin{:});
end

xlabel('f / Hz')


end