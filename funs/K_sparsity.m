function [ T ] = K_sparsity(x, K_s)
% This function find the threshold for K-Sparsity strategy

%% Input %%%%%%%%%%
%   x      : the input sequence
%   K_s    : the k-sparsity parameter (percentage) belongs to(0, 1)(default: 0.002)
%% Output %%%%%%%%%%
%   T      : the generated threshold
% Author : Zhibin Zhao
% Place  : Xi'an Jiaotong University
% Email  : zhibinzhao1993@gmail.com
% Date   : 2018.6
if nargin < 2
    K_s = 0.002;
end
k = round(length(x)*K_s);
Descend = sort(abs(x), 'descend');
T = Descend(k);

end

