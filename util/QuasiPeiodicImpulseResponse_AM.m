function Linchao  = QuasiPeiodicImpulseResponse_AM(N,fs)

t = (1:N)/fs;   


fc = 200;
fm = 10;
A1 = 0.0;
% Harmonic signal
y11 = A1; % 0.05*cos(2*pi*fm*t)+
y12 = cos(2*pi*fc*t);
y1  = y11.*y12;

% Periodic impulse response
f0 = 2000;             %调制频率
T0 = 0.002;           %延时
T  = 0.01;            %调制信号周期
AA  = 1;               %有效信号幅值

while (T0-T >= 0)
    T0 = T0-T;
end

zeta0R = 0.005;          %阻尼系数
zeta0L = 0.02;          %阻尼系数
sig11 = AA*exp(-zeta0R*(2*pi*f0*(t-T0)).^2).*cos(2*pi*f0*(t-T0));
sig11(1:(round(T0*fs)-1)) = 0;
sig12 = AA*exp(-zeta0L*(2*pi*f0*(t-T0)).^2).*cos(2*pi*f0*(t-T0));
sig12(round(T0*fs):N) = 0;
% sig11 = AA*exp(-(zeta0R/sqrt(1-zeta0R^2))*(2*pi*f0*(t-T0)).^2).*cos(2*pi*f0*(t-T0));
% sig11(1:(round(T0*fs)-1)) = 0;
%  sig12 = AA*exp(-(zeta0L/sqrt(1-zeta0L^2))*(2*pi*f0*(t-T0)).^2).*cos(2*pi*f0*(t-T0));
% sig12(round(T0*fs):N) = 0;
sig1 = sig11 + sig12;

% zeta0R = 0.1;          %阻尼系数
% zeta0L = 0.5;          %阻尼系数
% sig11 = AA*exp(-(zeta0R*(2*pi*f0*(t-T0)))/sqrt(1-zeta0R^2)).*sin(2*pi*f0*(t-T0));
% sig11(1:(round(T0*fs)-1)) = 0;
% sig12 = AA*exp((zeta0L*(2*pi*f0*(t-T0)))/sqrt(1-zeta0L^2)).*sin(2*pi*f0*(t-T0));
% sig12(round(T0*fs):N) = 0;
% sig1 = sig11 + sig12;

% zeta0R = 0.05;          %阻尼系数
% sig1 = AA*exp(-zeta0R*(2*pi*f0*(t-T0))) ./sqrt(1-zeta0R^2) .*cos(2*pi*f0*(t-T0));
% sig1(1:(round(T0*fs)-1)) = 0;

signal = sig1;
for irow = 1:N/(T*fs)
%     ImpulseIndex = round((2*rand(1)-1) * 0.1 * T *fs);
    ImpulseIndex = round((2*rand(1)-1) * 0.05 * T *fs);    
    Wss = round(T*fs*irow) + ImpulseIndex;
    for k = 1:N      
        if k > Wss            
            signal(k) = signal(k) + sig1(k-Wss);
        end
    end
end

y2 = signal;

Linchao = y1 + y2;
fm = 100/15;
% Linchao = Linchao * 0.5 .* (1-0.8*cos(2*pi*fm*t));