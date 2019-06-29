%% Demo1: 
% Author : Zhibin Zhao
% Place  : Xi'an Jiaotong University
% Email  : zhibinzhao1993@gmail.com
% Date   : 2018.6
clc
clear all
close all
addpath(genpath(fileparts(mfilename('fullpath'))));

Params = Config();
rng('default')
rng(Params.random_seed)
[ Sig , t] = Generate_Simulation(Params);


%% Perform the proposed method
[y_HHLP, cost_HHLP] = HHLP(Sig, Params);
%% Perform the square envelope spectrum kurtosis weighted L1-regularization method
Params.W_type = 'SESK';
Params.shrinkage = 'soft';
[y_L1_Kurtosis, cost_L1_Kurtosis] = HHLP(Sig, Params);
%% Perform the L1-regularization method
Params.W_type = 'None';
Params.shrinkage = 'soft';
[y_L1, cost_L1] = HHLP(Sig, Params);
%% Plot the results
figure(1)
subplot(411)
plot(t, Sig)
title('Original Signal')
ylabel('Amp (g)')
subplot(412)
plot(t, y_HHLP)
title('HHLP')
ylabel('Amp (g)')
subplot(413)
plot(t, y_L1_Kurtosis)
title('L1-Kurtosis')
ylabel('Amp (g)')
subplot(414)
plot(t, y_L1)
title('L1')
ylabel('Amp (g)')
xlabel('Time (s)')
filename = ['results', filesep, sprintf('Demo1_Extracted_Results.pdf')];
print(filename, '-dpdf');
%% Plot the history of the cost function
figure(2)
subplot(311)
plot(cost_HHLP)
title('HHLP')
subplot(312)
plot(cost_L1_Kurtosis)
title('L1-Kurtosis')
subplot(313)
plot(cost_L1)
title('L1')
xlabel('Iteration')
filename = ['results', filesep, sprintf('Demo1_Cost_Functions.pdf')];
print(filename, '-dpdf');
% If you feel our HHLP is useful for your research, please consider citing our paper:
% @article{zhao2019hierarchical,
%   title={Hierarchical hyper-Laplacian prior for weak fault feature enhancement},
%   author={Zhao, Zhibin and Wang, Shibin and An, Botao and Guo, Yanjie and Chen, Xuefeng},
%   journal={ISA transactions},
%   year={2019},
%   publisher={Elsevier}
% }