%% 
%gera��o do sinal
clear all
dt=0.001;
MC=100;
t=-500:dt:500;
x=zeros(1,length(t));
tau=-200:dt:200;
Rx=zeros(1,length(tau));
P=0;


for i=1:MC
    x=rtn_simple(t,dt,[0.1 1],1);
    avg=0;
    for j=1:length(t)
        avg=dt*x(j)+avg;
    end
    x=x-avg/1000;
    pin=0;
    for j=1:length(t)
        pin=dt*x(j)^2+pin;
    end
    P=P+pin/1000;
    for j=1:length(Rx)
        Rx(j)=x((length(x)+1)/2)*x(((length(x)+1)/2)-((length(tau)+1)/2)+j)+Rx(j);
    end
end

P/MC

%% 

SSD=fft(Rx/MC)/numel(tau);
FTsiga = double(abs(SSD)*2);    %Magnitude,Convert to double
%S_smth = sgolayfilt(FTsiga,20,501);

Fs = 1/dt;                          %Sampling frequency
Fn = Fs/2;                          %Nyquist frequency
Fv = linspace(0,1, fix(numel(tau)/2)+1)*Fn;  %Frequency vector

loglog(Fv,abs(SSD(1:length(Fv))))

sum=0;
for k=1:length(FTsiga)
sum=abs(PSD2(k))*0.0025+sum;
end
sum;
