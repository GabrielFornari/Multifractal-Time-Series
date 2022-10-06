function [beta, f, s] = psd(x, fflag)
% Power Spectrum Density of signal 'x'
% 
% Usage:
%   [beta] = psd(x, 1);
%
% Inputs:
%   x        signal (time series)
%   fflag    1|0 - for output log-log plot of s(pow) against f (freq)
%
% Outputs:
%   beta    Least Squares fit of log(s) against log(s)
%   f       Frequency values
%   s       Spectral density values
%

    if nargin < 2 
       fflag = 0;
    end
    
    n = length(x);

    % Ensure data has a mean of zero
    x = x - mean(x);

    % Welch window
    w = 1 - (([1:n]' - n/2)/(n/2)).^2;
    Wss = mean(w.^2);

    % calculate the spectrum
    k = n/2;
    s0 = (1/Wss)*(2/n)*abs(fft(w.*x)).^2;
    s = s0(1:k);
    f = [1:k]'/n;

    pol = polyfit(log(f), log(s), 1);
    beta = pol(1);
    
    if fflag == 1
        figure;
        loglog(f,s,'k.-');
        hold on;
        loglog(f, exp(polyval(pol, log(f))), 'k-', 'LineWidth', 2);
        hold off;
        xlabel('frequency (Hz)','FontSize',14);
        ylabel('PSD[|sfu|^2/Hz]','FontSize',14);
        title(['S(f) ~ f^{-\alpha} with \alpha = ' num2str(-pol(1)) ]);
    end
end
%  
%  Created by 
%       Gabriel Fornari
%  On 
%       20/05/2015 (dd/mm/yyyy)
%  
%  Based on
%       Copyright (c) 2005 Patrick E. McSharry (patrick@mcsharry.net)
%**********************************************************
