function [x] = mffiltering(H, sx, n, fflag)
% Modified Fourier filtering method to generate correlated series 
% with different behaviours.
% 
% Usage:
%   [x] = ff(H, sx, n, fflag);
%
% Inputs:
%   H        Hurst exponent (must be between 0.5 and 1.5)
%   sx       Scale crossover (fx = 1/sx)
%   n        length of the output series x
%   fflag    optional; 1|0 - for output plot of x
%
% Outputs:
%   x        signal (time series)
%

    if H < .5 || H > 1.5
        error('Parameter H must be between 0.5 and 1.5');
    end

    if nargin < 4
        fflag = 0;
    end
    
    beta = 2*H-1;
    
    n = 4*n;
    x = wgn(n, 1, 1);
    
    c = (fft(x));
    c(1) = [];
    m = 1:(length(c)/2)+1;
    c = c(m);
    f = m./n;
    
    fx = 1/sx;
    ix = find(f>fx);
    
    c(ix) = 2*c(ix)'.*((f(ix)/fx).^(-beta/2));
    
    ic = ifft(c);
    x = real(ic((n/8)+1:(n*3/8)));
    
    x = x-mean(x);
    
    if fflag
        figure;
        plot(x);
        title(strcat('Correlated Series (H=', num2str(H), ')'), ...
                     'FontSize', 24);
        xlabel('Time', 'FontSize', 20);
        ylabel('Amplitude', 'FontSize', 20);
        set(gca, 'FontSize', 16);
    end
end
