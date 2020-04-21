%% 
clc
clear all
close all
addpath(genpath(fileparts(mfilename('fullpath'))));
%% Figure initialization
FontSize = 11;   FontName = 'Times New Roman';
MarkerSize = 4;  LineWidth = 1;
%%
FlagFigureAutoSave = 1;
currentFolder = pwd;
rng('default')
rng(35)
%% Simulation
Fs = 20480;
N = 32768;
t = (0 : N-1) / Fs;
t = t(:);
% Sig_Cos = 1*cos(2*pi*200*t) + 1*cos(2*pi*4000*t);  lambda2 = 1.9; 
% harmonic interferences
Sig_Cos = (1 + 0.5*cos(2*pi*60*t)) .* cos(2*pi*700*t + 0.5*cos(2*pi*35*t))...
    + (1 + 0.5*cos(2*pi*120*t)) .* cos(2*pi*1400*t + 0.5*cos(2*pi*70*t)) + 1*cos(2*pi*150*t) + 1*cos(2*pi*550*t);


%% Setting the TQWT Parameters
Q = 2;
r = 5;
J = 10;
AH = @(Sig) tqwt_radix2(Sig, Q, r, J);
A = @(w) itqwt_radix2(w, Q, r , N);

%%
x_zero = AH(zeros(N,1));
RR = 100;
Error1 = zeros(RR, J + 1);
Error2 = zeros(RR, J + 1);
Error3 = zeros(RR, J + 1);
for i = 1 : RR
    fprintf(['Iteration:' num2str(i) '\n'])
    rng(i)
    Sig_Impulse = 2 * QuasiPeiodicImpulseResponse_AM(N,Fs);
    Sig_Impulse = Sig_Impulse(:);
    Noise = 1.2 * randn(N, 1);
    Sig_Combine = Sig_Cos + Sig_Impulse + Noise;
    x1 = AH(Sig_Combine);
    x2 = AH(Sig_Impulse);
    x3 = AH(Sig_Cos);
    for j = 1: J+1
        Error1(i, j) = sum(abs(x1{j}).^2);
        Error2(i, j) = sum(abs(x2{j}).^2);
        Error3(i, j) = sum(abs(x3{j}).^2);
    end
end


Error1 = bsxfun(@rdivide, Error1, max(Error1, [], 2));
Error2 = bsxfun(@rdivide, Error2, max(Error2, [], 2));
Error3 = bsxfun(@rdivide, Error3, max(Error3, [], 2));

% Directly calculate the Combined
Combined_ratio1 = mean(Error1, 1);
Combined_ratio_std1 = std(Error1, 1);

% Directly calculate the Impulse
Impulse_ratio1 = mean(Error2, 1);
Impulse_ratio_std1 = std(Error2, 1);

% Directly calculate the Harmonic
Cos_ratio1 = mean(Error3, 1);
Cos_ratio_std1 = std(Error3, 1);



%%
figure();clf;
ph(1) = bar(1:J+1, Impulse_ratio1, 'b');
% ph(2) = errorbar(1:J+1, Impulse_ratio1, Impulse_ratio_std1, 'r','LineWidth', LineWidth+1, 'Marker', 'none', 'LineStyle', 'none');
box off
xlim_min = 0; xlim_max = J+2;
ylim_min = 0;            ylim_max = 1.2;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);

filename = ['Results', filesep, sprintf('Impulse_Distribution.pdf')];
print(filename, '-dpdf');

%%
figure();clf;
ph(1) = bar(1:J+1, Cos_ratio1, 'b');
% ph(2) = errorbar(1:J+1, Cos_ratio1, Cos_ratio_std1, 'r','LineWidth', LineWidth+1, 'Marker', 'none', 'LineStyle', 'none');
box off
xlim_min = 0; xlim_max = J+2;
ylim_min = 0;            ylim_max = 1.2;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);




filename = ['Results', filesep, sprintf('Harmonic_Distribution.pdf')];
print(filename, '-dpdf');