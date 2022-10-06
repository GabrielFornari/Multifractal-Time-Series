function [alpha, falpha, Dq, tauq] = abmp(p, q, fflag)
% Calculates the multifractal spectrum according to binomial
% multiplicative process.
% 
% Usage:
%   [alpha, falpha, Dq, tauq] = abmp(p, q, fflag)
%
% Inputs:
%   p        eddy partition measure (0 <= p <= 1)
%   q        weight exponent
%   fflag    optional; 1|0 - for output plot of falpha versus alpha
%
% Outputs:
%   alpha    singularity strength
%   falpha   Hausdorff dimension
%   Dq       Fractal dimension
%   tauq     mass exponent
%
    
    if nargin < 3
        fflag = 0;
    end

    idx = find(q==1);
    
    tauq = -log(p.^q + ((1-p).^q))/log(2);
    Dq = -tauq./(1-q);
    
    if ~isempty(q)
        qtau = q;
        qtau(idx) = [];
        Dq(idx) = [];
    end
    
    % First Derivative - Legendre Transform
    alpha = diff(tauq);
    for iAlpha = 1:length(tauq)-1
        alpha(iAlpha) = alpha(iAlpha)/(q(iAlpha+1)-q(iAlpha));
    end
    falpha = (q(1:end-1).*alpha) - tauq(1:end-1);
    
    if fflag
        figure;
        plot(qtau, Dq, 'LineWidth', 2.5);
        title('D_{q} versus q', 'FontSize', 24);
        xlabel('q', 'FontSize', 20);
        ylabel('Dq', 'FontSize', 20);
        set(gca, 'FontSize', 16);
        
        figure;
        plot(qtau, alpha, 'LineWidth', 2.5);
        title('Lipschitz Holder exponent versus q', 'FontSize', 24);
        xlabel('q', 'FontSize', 20);
        ylabel('\alpha', 'FontSize', 20);
        set(gca, 'FontSize', 16);
        
        figure
        plot(q, tauq, 'LineWidth', 2.5);
        title('\tau(q) versus q', 'FontSize', 24);
        xlabel('q', 'FontSize', 20);
        ylabel('\tau(q)', 'FontSize', 20);
        set(gca, 'FontSize', 16);
        
        figure;
        plot(alpha, falpha, 'LineWidth', 2.5)
        title('Singularity Spectrum - f(\alpha) versus \alpha', ...
                               'FontSize', 24);
        xlabel('\alpha', 'FontSize', 20);
        ylabel('f(\alpha)', 'FontSize', 20);
        set(gca, 'FontSize', 16);
    end
end
