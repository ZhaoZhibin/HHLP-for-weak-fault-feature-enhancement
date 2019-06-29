function [ Params ] = Config()
% This function performs parameter configuration
% Author : Zhibin Zhao
% Place  : Xi'an Jiaotong University
% Email  : zhibinzhao1993@gmail.com
% Date   : 2018.6

%% Set the random seed to make sure the reproducibility
Params.random_seed = 36;          % The random state

%% Parameters of Generating Simulation
Params.Fs            = 20480;     % The sampling frequency of the simulation signal
Params.N             = 4096;      % The length of the signal
Params.mixture_ratio = [1, 1.2, 1.2]; % The mixing ratio of [impulses, harmonic, noise].
% Periodic Impulses
Params.AP         = 2;            % The amplitude of the impulse
Params.T          = 0.01;         % The fault period
Params.tau0       = 0.002;        % The initial phase
Params.f1         = 2000;         % The resonance frequency
Params.zeta_left  = 0.02;         % The left damping
Params.zeta_right = 0.005;        % The right damping
% Harmonic interference
Params.Order = 2;                 % The Order of the interference
Params.CF    = 700;               % The basic carrier frequency
Params.AM    = 60;                % The amplitude modulation frequency 
Params.FM    = 35;                % The frequency modulation frequency
Params.H     = [150, 550];        % The discrete frequencies
% noise type
Params.noise_type = 'Gaussian';   % The noise type can be 'Gaussian' or 'Laplacian'


%% Parameters of the HHLP
% TQWT parameters
Params.Q         = 2;            % The Q factor(suggesting a small value 1<Q<5; default: 2)
Params.r         = 5;            % The redundant factor (default: 5)
Params.J         = 10;           % The level factor (default: 10)
Params.algo_type = 'radix2';     % The algorithm of TQWT: 'radix2' or 'None' (default: 'radix2')
% Weight parameters
Params.W_type     = 'multi-scale PMI';      % type:'SESK', 'multi-scale PMI' or 'None'
Params.F_Interval = [96, 104];              % The searching interval of the characteristic frequency
% Optimization parameters
Params.K_s       = 0.002;        % The K-sparsity parameter (percentage) belongs to(0, 1) (default: 0.002)
Params.shrinkage = 'half';       % The shrinkage function: 'soft' or 'half'(default: 'half')
Params.Nit       = 100;          % Number of iterations (default: 100)

end

