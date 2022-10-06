function [f, P] = powspec(varargin)
% powspec    Power Spectrum of a series X
%   This function estimates the frequency and the power of a signal 
%   through the Fast Fourier Transform (FFT). The spectrum is 
%   calculated by squaring the Fourier coefficients.
%
% Usage:
%   [f, P] = powspec(t, x, , 'velocity', v, window, plot)  
%                                      calculates the power(P) for 
%                                      each frequency(F) of a time 
%                                      series(X)
%
% Inputs:
%   t      time (vector with the time of each sample)
%   x      time series (vector with samples)
%   v      optional; mean velocity during the measure, must be used 
%          after 'velocity'
%   win    optional; can be:
%             'hanning'  - multiplies the time series by a hanning 
%                          window
%             'hamming'  - multiplies the time series by a hamming 
%                          window
%             'blackman' - multiplies the time series by a blackman 
%                          window
%   'plot' optional; plot the spectrum in a new figure
%
% Outputs:
%   f    frequency
%   P    power
%

    if nargin < 2
        error('Insuficient number of parameters');
    else
        t = varargin{1};
        t = t(:);
        x = varargin{2};
        x = x(:);
        
        N = length(x);
        dt = (t(end)-t(1))/(N-1);
        fs = 1/dt;
    end
    if any(strcmpi(varargin, 'hanning'))
        win = hanning(N);
        x=x.*win;
    elseif any(strcmpi(varargin, 'hamming'))
        win = hamming(N);
        x=x.*win;
    elseif any(strcmpi(varargin, 'blackman'))
        win = blackman(N);
        x=x.*win;
    else
        win = ones(N, 1);
    end

    S = sum(win);
    c = fft(x);
    c = c(2:(N/2+1));
    P = (2*abs(c).^2)./(S^2);

    f = (1:N/2)*(fs/N);

    xLabelName = 'Frequency (Hz)';

    if any(strcmpi(varargin, 'velocity'))
        iArg = 1;
        while ~strcmpi(varargin(iArg), 'velocity')
            iArg = iArg+1;
        end
        vel = varargin{iArg+1};
        f = f.*(2*pi)/vel;
        xLabelName = 'Wavenumber (rad/km)';
    end

    if any(strcmpi(varargin, 'plot'))
        figure;
        subplot(2,1,1);
        plot(t, x,'k-');
        title('Time Series');
        subplot(2,1,2);
        loglog(f, P, 'k-');
        title('Power Spectrum');
        xlabel(xLabelName);
    end
end
