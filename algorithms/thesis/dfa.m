function [alpha, f] = dfa(x, scales, order, fflag)
% Detrended Fluctuation Analysis of signal 'x'.
% 
% Usage:
%   [alpha, f] = dfa(x, scales, order, fflag);
%
% Inputs:
%   x        signal (time series)
%   scales   vector with all desired scales
%   order    optional; polynomial order for the least square 
%            (order = 1 standard DFA)
%   fflag    optional; 1|0 - for output log-log plot of F('scales')
%                            versus 'scales'
%
% Outputs:
%   alpha    fit of log(f) versus log(scales)
%   f        average fluctuation
%

    if nargin < 4
        fflag = 0;
    end
    if nargin < 3
        order = 1;
    end
    
    xLen = length(x);
    
    %Step 1 - Cumulative sum minus mean
    y = cumsum(x(:) - mean(x));
    %y = (x(:) - mean(x));
    
    nScales = length(scales);
    f = zeros(nScales, 1);
    for iScale=1:nScales
        scale = scales(iScale);
        
        %Step 2 - Divide 'y' into 'n' segments of length 's'
        nSgm = floor(xLen/scale);

        %Step 3 - Calculate local trend for each segment by 
        %                                          least-square
        fq = zeros(nSgm, 1);
        for iSgm=1:nSgm
            values = scale*(iSgm-1)+1:scale*iSgm;
            pol = polyfit(values, y(values)', order);
            fq(iSgm) = (1/scale)*sum((y(values)' - ...
                                       polyval(pol, values)).^2);
        end
        fqInv = zeros(nSgm, 1);
        for iSgm=1:nSgm
            values = xLen-scale*iSgm+1:xLen-scale*(iSgm-1);
            pol = polyfit(values, y(values)', order);
            fqInv(iSgm) = (1/scale)*sum((y(values)' - ...
                                       polyval(pol, values)).^2);
        end    
        f(iScale) = sqrt((sum(fq)+sum(fqInv))/(2*nSgm));
    end
    
    pol = polyfit(log(scales)', log(f), 1);
    alpha = pol(1);
    
    if fflag
        figure;
        loglog(scales, f, 'ko', 'MarkerSize', 10);
        hold on;
        loglog(scales, exp(polyval(pol, log(scales))), ...
                                        'k-', 'LineWidth', 2);
        title(['F_{2}(s) with h_{2} = ' num2str(alpha)], ...
                                        'FontSize', 24);
        xlabel('scales', 'FontSize', 20);
        ylabel('F_{2}(s)', 'FontSize', 20);
        set(gca, 'FontSize', 16);
    end
end
