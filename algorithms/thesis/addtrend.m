function [px] = addtrend(x, a, p, fflag)
% Combine a monotonous trend with a signal by adding polynomials of
% different order p to the time series x.
% 
% Usage:
%   [px] = ff(x, a, p, fflag);
%
% Inputs:
%   x        time series
%   a        a constant (e.g. 10^2)
%   p        polynomial order (not necessarily an integer)
%   fflag    optional; 1|0 - for output plot of px
%
% Outputs:
%   px       modified signal (time series)
%

    x = x(:);
    n = length(x);
    i = (1:n)'/n;
    
    length(x)
    length(i)
    
    px = x + a.*(i.^p);
    
    if fflag
        figure;
        plot(px);
        title(strcat('Series with trend (a=', num2str(a), ', p=', ...
                     num2str(p), ')'), 'FontSize', 24);
        xlabel('Time', 'FontSize', 20);
        ylabel('Amplitude', 'FontSize', 20);
        set(gca, 'FontSize', 16);
    end
end
