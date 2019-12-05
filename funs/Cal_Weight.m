function [ Weight ] = Cal_Weight(x, A, AH, Params)
% This function calculate the weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x:                          The input signal
% A:                          Function handles for the transform A
% AH:                         Function handles for the inverse transform A
% Params: a struct contains all the parameters
%       Params.W_type:        The weight type: 'SESK', 'multi-scale PMI' or 'None'
%                             SESK: square envelope spectrum kurtosis (default)
%                             multi-scale PMI: multi-scale periodic modulation intensity
%                             None: the weight is inoperative
%       Params.Fs:            The sampling frequency of the simulation signal
%       Params.F_Interval:    The searching interval of the characteristic frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weight:                     The generated weight vector

% Author : Zhibin Zhao
% Place  : Xi'an Jiaotong University
% Email  : zhibinzhao1993@gmail.com
% Date   : 2018.6

if ~isfield(Params, 'W_type')
    W_type = 'SESK';
else
    W_type = Params.W_type;
end
fprintf(sprintf('HHLP with W_type: %s.\n', W_type));
% Define the weight information
J1 = numel(x);
N = length(A(x));
switch W_type
    case 'None'
        Weight = ones(J1, 1);
    case 'multi-scale PMI' 
        if ~isfield(Params, 'Fs')
            error('Multi-scale PMI need the sampling frequency');
        end
        if ~isfield(Params, 'F_Interval')
            error('Multi-scale PMI need the searching interval of characteristic frequency');
        end
        x_tmp  = x;
        K = zeros(J1, 1);
        for i = 1: J1
            Hx = abs(hilbert(x_tmp{i}));
            Ratio = N / length(x_tmp{i});
            fs = Params.Fs / Ratio;
            [PMI, ~, ~] = myHNR (Hx, fs, Params.F_Interval);
            K(i) = abs(PMI);
        end     
        Weight = max(K) ./ K;   
    case 'SESK'
        if ~isfield(Params, 'Fs')
            error('The SESK method need the sampling frequency');
        end       
        x_zero = AH(zeros(N,1));
        x_tmp  = x;
        Kurtosis = zeros(J1, 1);
        for i = 1: J1
            Tmp3 = x_zero;
            Tmp3{i} = x_tmp{i};
            y3 = real(A(Tmp3));
            [ yf3, ~ ] = Hilbert_envelope( y3 , Params.Fs , 1);
            Kurtosis(i) = kurtosis(yf3);      
        end 
        Weight = max(Kurtosis) ./ Kurtosis;       
    otherwise
        error('Unknown method.')
end      
end

