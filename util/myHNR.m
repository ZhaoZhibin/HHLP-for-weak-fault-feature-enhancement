function [HNR,f0,corr_curve] = myHNR (x,fs,rg,opt)
% This function is directly from Ming Zhao
% *************************************************************************
% Calculate the HNR of a real signal
% Ming Zhao @ UC, April.13, 2016
% ---------------------------INPUT-----------------------------------------
% x -------- Input signal;
% fs ------- Sampling frequency in Hz,
% rg ------- frequency search range in Hz,default: [0 fs/2]
% opt ------ specify the which method you use,default:'unbiased'
% ---------------------------OUTPUT----------------------------------------
% HNR ------ The returned HNR value
% f0  ------ Fundamental frequency corresponding to the detected period
% corr_curve ----- The correlation curve
% *************************************************************************

% Substract the mean from the signal
y = x - mean(x);

% Calculate the correlation curve
% NA = xcorr(y);
if nargin <= 3
    NA = xcorr(y,'unbiased');
else
    NA = xcorr(y,opt);
end
  
%
% Since NA is symmetric, only the Last half of correlation curve is
% extracted for further analysis
NA = NA( ceil(length(NA)/2):end );

% find first zero crossing
sample1 = NA(1);

for lag = 2:length(NA)
    sample2 = NA(lag);
    if(( sample1 > 0 ) && ( sample2 < 0 ))
        zero_position = lag;
        break;
    elseif(( sample1 == 0 ) || ( sample2 == 0))
        zero_position = lag;
        break;
    else
        sample1 = sample2;
    end
end

% Return the correlation curve
corr_curve = NA;

% Discard the samples before the first zero-crossing point
ZA = NA( zero_position:end );

%  Find the maximum value and its position
[max_value max_position] = max(ZA);
max_position = max_position+zero_position-1;

if nargin == 3
    ind_max = round(fs/rg(1)+1);
    ind_min = round(fs/rg(2)+1);
    % Discard the samples before the first zero-crossing point
    ZA = NA( ind_min:ind_max );
    
    %  Find the maximum value and its position
    [max_value max_position] = max(ZA);
    max_position = max_position+ind_min-1;
 
end

% Return the Fundamental frequency
f0 = fs/(max_position-1);

% Returned HNR value
HR = max_value;
HNR= HR/(NA(1)-HR);

end
    
    
    
    
    
     
    
    
    
