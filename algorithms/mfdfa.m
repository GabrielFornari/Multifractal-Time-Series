function [alpha, falpha, q, tauq] = mfdfa(x, scales, q, order, fflag)
% Calculates the multifractal spectrum of signal 'x' using DFA method
% 
% Usage:
%   [alpha, falpha, q, tauq] = mfdfa(x, [4:64], [-20:20], 1, 1)
%
% Inputs:
%   x        signal (time series)
%   scales   vector with all desired scales
%   q        vector with all q-order that weight the variations
%   order    polynomial order for the least square (order = 1 standard DFA)
%   fflag    1|0 - for output plot of singularity spectrum and vector of
%                  moments
%
% Outputs:
%   alpha    Singularity Strength
%   falpha   Hausdorff Dimension
%   q        List of Exponents
%   tauq     Vector of Moments
%
    if nargin < 5
        fflag = 0;
    end
    if nargin < 4
        order = 1;
    end
    
    nq = length(q);
    xLen = length(x);
    
    %Step 1 - Cumulative sum minus mean
    y = cumsum(x(:) - mean(x));
    
    nScales = length(scales);
    f = zeros(nq, nScales);
    for iScale=1:nScales
        scale = scales(iScale);
        
        %Step 2 - Divide 'y' into 'n' segments of length 's'
        nSgm = floor(xLen/scale);

        %Step 3 - Calculate local trend for each segment by least-square
        fq = zeros(nSgm, 1);
        for iSgm=1:nSgm
            values = scale*(iSgm-1)+1:scale*iSgm;
            pol = polyfit(values, y(values)', order);
            fq(iSgm) = (1/scale)*sum((y(values)' - polyval(pol, values)).^2);
        end
        fqInv = zeros(nSgm, 1);
        for iSgm=1:nSgm
            values = xLen-scale*iSgm+1:xLen-scale*(iSgm-1);
            pol = polyfit(values, y(values)', order);
            fqInv(iSgm) = (1/scale)*sum((y(values)' - polyval(pol, values)).^2);
        end    
        for iq = 1:nq
            if(q(iq) == 0)
                f(iq, iScale) = exp(1/(4*nSgm) * (sum(log(fq))+sum(log(fqInv))));
            else
                f(iq, iScale) = (1/(2*nSgm)*(sum(fq.^(q(iq)/2))+sum(fqInv.^(q(iq)/2))))^(1/q(iq));
            end
        end
    end
    h = zeros(nq, 1);
    for iq = 1:nq
        coef = polyfit(log(scales), log(f(iq, :)), 1);
        h(iq) = coef(1);
    end
    
    tauq = h.*q' - 1;
    % First Derivative - Legendre Transform
    alpha = diff(tauq);
    for iAlpha = 1:length(tauq)-1
        alpha(iAlpha) = alpha(iAlpha)/(q(iAlpha+1)-q(iAlpha));
    end
    falpha = (q(1:end-1).*alpha') - tauq(1:end-1)';
    
    % Plots
    if fflag
        figure;
        plot(q, tauq, 'ko', 'MarkerSize', 10);
        title('\tau(q) versus q', 'FontSize', 24);
        xlabel('q', 'FontSize', 20);
        ylabel('\tau(q)', 'FontSize', 20);
        set(gca, 'FontSize', 16);
        
        figure;
        plot(alpha, falpha, 'kx', 'MarkerSize', 10);
        title('Singularity Spectrum - f(\alpha) versus \alpha', 'FontSize', 24);
        xlabel('\alpha', 'FontSize', 20);
        ylabel('f(\alpha)', 'FontSize', 20);
        set(gca, 'FontSize', 16);
    end
end
