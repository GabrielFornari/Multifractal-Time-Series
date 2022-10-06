function [ts, wd] = randStaircase(steps, weights, seed, fflag)
% Creates a random discrete staircase (see Muzy, 1994)
%
% Usage
%   [ts, wd] = randStaircase(7, [.5, 0, .5], 0, 1);
%
% Inputs
%   steps    number of loops
%   weights  vector with weights
%              [.5, 0, .5] -> Monofractal Staircase
%              [.69, .46, -.46, .31] -> Multifractal Staircase
%   seed     seed for random number generator
%   fflag      1|0 - for output plot
%
% Outputs
%   ts       discrete time series (length ~ size(weights)^steps)
%   wd       weight distribution (length ~ size(weights)^steps)
%
    if nargin < 4
        fflag = 0;
    end
    if nargin > 2
        rng(seed);
    end
    
    base = length(weights);
    n = base^steps;
    
    if mod(base, 2) % if 'base' is odd
        n = n-1;
    end
    
	wd = ones(n, 1);
	time = (0:(n-1))./n;
	iTime = 0;

    while(iTime < steps)
        divs = floor(base.*time);
        time   = base .*time - divs;
       
        weights = weights(randperm(base));
        
        wd = wd.*weights(divs+1)';
        iTime = iTime+1;
    end
    
	ts = cumsum(wd);
    
    if fflag == 1
        figure;
        plot(ts, 'k');
        title('Generalized Random Staircase', 'FontSize', 24);
        xlabel('x', 'FontSize', 20);
        ylabel('f(x)', 'FontSize', 20);
        set(gca, 'FontSize', 16);
    end
end

%  
%  Created by 
%       Gabriel Fornari (gabriel.fornari@inpe.br)  
%  At 
%       18/05/2015 (dd/mm/yyyy)
%  
%  Based on
%       Part of Wavelab Version 850
%       Built Tue Jan  3 13:20:39 EST 2006
%       This is Copyrighted Material
%       For Copying permissions see COPYING.m
%       Comments? e-mail wavelab@stat.stanford.edu
