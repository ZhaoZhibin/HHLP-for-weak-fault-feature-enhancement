function [y_out, cost] = HHLP(y, Params)
% This function solves the following optimization problem:
% with (p=1: soft thresholding, p=0.5: half thresholding)
% min 1/2 * ||y - A*x||^2 + lambda * ||Wx||_p^p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y:                          The input signal
% Params: a struct contains all the parameters
%%% TQWT parameters
%       Params.Q:             The Q factor(suggesting a small value 1<Q<5; default: 2)
%       Params.r:             The redundant factor (default: 5)
%       Params.J:             The level factor (default: 10)
%       Params.algo_type:     The algorithm of TQWT: 'radix2' or 'None' (default: 'radix2')
%%% Optimization parameters
%       Params.K_s:           The K-sparsity parameter (percentage) belongs to(0, 1) (default: 0.002)
%       Params.shrinkage:     The shrinkage function: 'soft' or 'half'(default: 'half')
%       Params.Nit:           Number of iterations (default: 100)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y_out:                      The denoised signal
% cost:                       The history of the cost function


% Reference: 'Hierarchical hyper-Laplacian prior for weak fault feature
% enhancement', ISA Transactions, 2019
% https://zhaozhibin.github.io/
% Author : Zhibin Zhao
% Place  : Xi'an Jiaotong University
% Email  : zhibinzhao1993@gmail.com
% Date   : 2018.6

N = length(y);
y_means = mean(y);
y = y - y_means;
if nargin < 2
    Params = struct([]);
end
% initialization
if ~isfield(Params, 'Q')
    Q = 2;
else
    Q = Params.Q;
end
if ~isfield(Params, 'r')
    r = 5;
else
    r = Params.r;
end
if ~isfield(Params, 'J')
    J = 10;
else
    J = Params.J;
end
if ~isfield(Params, 'algo_type')
    algo_type = 'radix2';
else
    algo_type = Params.algo_type;
end
if ~isfield(Params, 'K_s')
    K_s = 0.002;
else
    K_s = Params.K_s;    
end
if ~isfield(Params, 'shrinkage')
    shrinkage = 'half';
else
    shrinkage = Params.shrinkage;     
end
if ~isfield(Params, 'Nit')
    Nit = 100;
else
    Nit = Params.Nit;      
end
% Define the transformation and the initialization
if strcmp(algo_type, 'radix2')
    if rem(log2(N), 1)
        error('If we use radix2, the length of y should be the power of 2');
    end
    AH = @(Sig) tqwt_radix2(Sig, Q, r, J);
    A = @(w) itqwt_radix2(w, Q, r , N);
    normA = ComputeNow(N, Q, r, J, 'radix2');
else
    AH = @(Sig) tqwt(Sig, Q, r, J);
    A = @(w) itqwt(w, Q, r , N);
    normA = ComputeNow(N, Q, r, J);  
end
fprintf('\n');
%% starting HHLP
fprintf(sprintf('performing starting HHLP...\n'));
% Initialize the coefficients
init = AH(y);

% Define the weight information
Weight = Cal_Weight(init, A, AH, Params);
fprintf(sprintf('HHLP with shrinkage: %s.\n', shrinkage));
% Other initialization
cost  = zeros(Nit , 1);
mu    = 0.9;
x     = init;
x_old = init;
AHy   = AH(y);
Ax    = A(x);
iter  = 1;
% Define the operator
half  = @(x, T) (2/3*x.*( 1 + cos( 2*pi/3 - 2/3*acos(T/8*(abs(x)/3).^(-3/2)) ) )) .* (abs(x)>54^(1/3)/4*T^(2/3));
soft  = @(x, T) max(1 - T./abs(x), 0) .* x;
while iter <= Nit
%     fprintf(['Iteration:' num2str(iter) '\n'])
    % forward operator
    % x = x - mu * (A^T(A(x) - y))
    tmp = x;
    AHAx = AH(Ax);    
    for i = 1:numel(tmp)
        tmp{i} = tmp{i} + mu * ( AHy{i} - AHAx{i} );
    end   
    % backward (thresholding) with K-sparsity
    Temp = [];
    for i = 1: J+1   
        Temp = [Temp ; tmp{i}(:) / normA(i) / Weight(i)];
    end
    T = K_sparsity(Temp, K_s);
    for i = 1: numel(x)
        switch shrinkage
            case 'half'
                x{i} = half(tmp{i}, T * normA(i) * Weight(i));
            case 'soft'
                x{i} = soft(tmp{i}, T * normA(i) * Weight(i));
            otherwise
                error('Undefined shrinkage.')
        end
    end
    
    Ax = A(x);
    % cost function history
    cost(iter) = 0.5 * norm(y(:)-Ax(:))^2 ;  
    stop = 0;
    total = 0;
     for i = 1: J+1
        switch shrinkage
            case 'half'
                cost(iter) = cost(iter) + Weight(i) * sum(abs(x{i}(:)).^0.5);
            case 'soft'
                cost(iter) = cost(iter) + Weight(i) * sum(abs(x{i}(:)));
        end         
        stop = stop + sum((x{i} - x_old{i}).^2);
        total = total + sum((x{i}).^2);
    end   
    % checkpoint
    if sqrt(stop) / sqrt(total) < 1e-6
        break;
    end
    iter = iter + 1;
    x_old = x;
end
if iter < Nit
    cost = cost(1:iter);
end
y_out = real(A(x)) + y_means;
y_out = y_out(:);
fprintf(sprintf('HHLP finished...\n'));
end




