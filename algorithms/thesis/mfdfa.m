function [alpha, fAlpha, q, tauq, hq] = mfdfa(x, scales, q, order, fflag)
% Calculates the multifractal spectrum of signal 'x' using DFA method
% 
% Usage:
%   [alpha, falpha, q, tauq] = mfdfa(x, [4:64], [-20:20], 1, 1);
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
    hq = zeros(nq, 1);
    for iq = 1:nq
        coef = polyfit(log(scales), log(f(iq, :)), 1);
        hq(iq) = coef(1);
    end
    
    tauq = hq.*q' - 1;
     
    % First Derivative - Legendre Transform
    alpha = diff(tauq);
    for iAlpha = 1:length(tauq)-1
        alpha(iAlpha) = alpha(iAlpha)/(q(iAlpha+1)-q(iAlpha));
    end
    fAlpha = (q(1:end-1).*alpha') - tauq(1:end-1)';
    
    % Remove discontinuity
    q1 = q;
    idx = (q>.5 & q<1.5); %(q > 0 & q < 2);
    q1(idx) = [];
    
    Dq = tauq./(q-1)';
    Dq(idx) = [];
    
    % Plots
    if fflag      
        figure;
        loglog(scales, f, '.', 'MarkerSize', 20, 'LineWidth', 2.5);
        title('F_{q} versus scales', 'FontSize', 24);
        xlabel('scales', 'FontSize', 20);
        ylabel('F_{q}', 'FontSize', 20);
        legend(num2str(q'));
        set(gca, 'FontSize', 16);
        hold on;
        for iq = 1:nq
            loglog(scales, exp(polyval(polyfit(log(scales), log(f(iq, :)), 1), log(scales))));
        end
        hold off;
        
        figure;
        plot(q, hq, 'ko', 'MarkerSize', 6, 'LineWidth', 2.5);
        title('Hurst Exponent h(q) versus q', 'FontSize', 24);
        xlabel('q', 'FontSize', 20);
        ylabel('h(q)', 'FontSize', 20);
        set(gca, 'FontSize', 16);
        
        figure;
        plot(q1, Dq, 'ko', 'MarkerSize', 10, 'LineWidth', 2.5);
        title('Generalized Dimensions D(q) versus q', 'FontSize', 24);
        xlabel('q', 'FontSize', 20);
        ylabel('Dq', 'FontSize', 20);
        set(gca, 'FontSize', 16);
        
        figure;
        plot(q, tauq, 'ko', 'MarkerSize', 10, 'LineWidth', 2);
        title('Mass exponent - \tau(q) versus q', 'FontSize', 24);
        xlabel('q', 'FontSize', 20);
        ylabel('\tau(q)', 'FontSize', 20);
        set(gca, 'FontSize', 16);
        
        figure;
        plot(alpha, fAlpha, 'ko', 'MarkerSize', 10, 'LineWidth', 2);
        title('Singularity Spectrum - f(\alpha) versus \alpha', 'FontSize', 24);
        xlabel('\alpha', 'FontSize', 20);
        ylabel('f(\alpha)', 'FontSize', 20);
        set(gca, 'FontSize', 16);
    end
end
