function [alpha, falpha] = apmodel(p, l1, l2, fflag, n)
% Calculates the multifractal spectrum according to p-model.
% 
% Usage:
%   [alpha, falpha] = apmodel(p, l1, l2, fflag, n)
%
% Inputs:
%   p        eddy partition measure (0 <= p <= 1)
%   l1       eddy scale 1
%   l2       eddy scale 2
%   fflag    optional; 1|0 - for output plot of falpha versus alpha
%   n        optional; number of points of alpha
%
% Outputs:
%   alpha    singularity strength
%   falpha   Hausdorff dimension
%
    
    if nargin < 4
        fflag = 0;
    end
    if nargin < 5
        n = 100;
    end
    if p > 1
        error('Arg p have to be beetween 0 and 1');
    end
    
    m = 2:n-1;
    p1 = p;
    p2 = 1-p;
    
    falpha = zeros(n, 1);
    alpha = zeros(n, 1);
    
    falpha(end-1:-1:2) = ((n./m-1).*log(n./m-1)-(n./m).*...
                         log(n./m))./(log(l1)+(n./m-1).*log(l2));
    alpha(end-1:-1:2) = (log(p1)+(n./m-1)*log(p2))./...
                        (log(l1)+(n./m-1)*log(l2));
    
    % leftmost point
    alpha(1) = log(p1)/log(l1);
    
    % rightmost point
    alpha(end) = log(p2)/log(l2);
    
    if fflag
        figure;
        plot(alpha, falpha);
        title('Singularity Spectrum', ...
              'FontSize', 24);
        xlabel('\alpha', 'FontSize', 20);
        ylabel('f(\alpha)', 'FontSize', 20);
        set(gca, 'FontSize', 16);
    end
end
