function [scales] = logscales(a, b, n)
% Generates exponentially spaced numbers beetween lower and 
% upper limits.
%
% Usage:
%   [scales] = pmodelSignal(a, b, n);
%
% Inputs:
%   a       lower bound; first scale
%   b       upper bound; last scale
%   n       optional; number os scales (if not defined, 
%                                   n = (log2(b)-log2(a))*4)
%
% Outputs:
%   scales  array with exponentially spaced scales
%
    
    if nargin < 3
       n = (log2(b) - log2(a))*4;
    end

    scales = round(logspace(log10(a), log10(b), n));
    scales = unique(scales, 'first');
end
