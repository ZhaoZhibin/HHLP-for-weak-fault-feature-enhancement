%% 
clc
clear all
close all
addpath(genpath(fileparts(mfilename('fullpath'))));
FontSize = 11;   FontName = 'Times New Roman';
MarkerSize = 7;  LineWidth = 1;
%%
FlagFigureAutoSave = 1;
currentFolder = pwd;

%% Simulation
Fs = 20480;
N = 409600;
mode = 'outer';
t = (0 : N - 1) / Fs;
t = t(:);
Sig_Cos = (1 + 0.5*cos(2*pi*60*t)) .* cos(2*pi*700*t + 0.5*cos(2*pi*35*t))...
    + (1 + 0.5*cos(2*pi*120*t)) .* cos(2*pi*1400*t + 0.5*cos(2*pi*70*t)) + 1*cos(2*pi*150*t) + 1*cos(2*pi*550*t);


%% Setting the TQWT Parameters
Q = 2;
r = 5;
J = 10;
AH = @(Sig) tqwt_radix2(Sig, Q, r, J);   
A = @(w) itqwt_radix2(w, Q, r , N);
now = ComputeNow(4096,Q,r,J,'radix2');
Temp = [];
Energy = zeros(100, J+1);
%% 把信号分成100组计算小波系数，然后进行概率拟合
for j = 1 : 100
    x = AH(Sig_Cos((j-1)*4096+1:j*4096));
    Temp1 = [];
    for i = 1:numel(x)
        Temp1 = [Temp1; x{i}(:) / now(i)];
        Energy(j, i) = norm(x{i}(:) / now(i));
    end    
    Temp = [Temp; Temp1];
end
Temp = Temp - mean(Temp);


%% Fitting the probability distribution
[Number,edges] = histcounts(Temp, 2000, 'Normalization', 'probability');
for i = 1 : length(Number)
    Points(i) = (edges(i + 1) + edges(i))/2;
end

%% Calculate the probability
% plot(Points,log2(N),'r-','LineWidth',1.5)
% axis([-0.4, 0.4, -15, 0])
% hold on
y = -1.5:0.01:1.5;
g = max(Number);
% gaussian
k1 = 40;
f1 = exp(-abs(y).^2.*k1) * (g);
% laplacian
k2 = 18;
f2 = exp(-abs(y).^1.*k2) * (g);
% hyper-laplacian 0.2
k3 = 10;  
f3 = exp(-abs(y).^0.1.*k3) * (g);
% hyper-laplacian 0.5
k4 = 12;     %16
f4 = exp(-abs(y).^0.5.*k4) * (g);  

%% Print the time domain

figure();clf;
hold on
ph(1) = plot(Points, log2(Number), 'b-*', 'LineWidth', LineWidth + 1,'MarkerIndices',796);%
ph(2) = plot(y, log2(f1), 'k->', 'LineWidth', LineWidth + 1,'MarkerIndices',130);
ph(3) = plot(y, log2(f2), 'g-o', 'LineWidth', LineWidth + 1,'MarkerIndices',130);
ph(4) = plot(y, log2(f3), 'r-^', 'LineWidth', LineWidth + 1,'MarkerIndices',140);
ph(5) = plot(y, log2(f4), 'm-d', 'LineWidth', LineWidth + 1,'MarkerIndices',130);
hold off
box off
legend1 = legend(ph, 'Empirical' , 'Gaussian(p=2)', 'Laplacian(p=1)', 'Hyper-Laplacian(p=0.2)', 'Hyper-Laplacian(p=0.5)');
set(legend1,'position',[0.361909575703485 0.15911596451407 0.528089875334435 0.329597259975321],'Orientation','vertical', 'FontSize',FontSize,'FontName',FontName)
legend boxoff

xlim_min = -0.6; xlim_max = 0.6;
ylim_min = -20;            ylim_max = 0;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);

filename = ['Results', filesep, sprintf('Harmonic_Probability.pdf')];
print(filename, '-dpdf');
