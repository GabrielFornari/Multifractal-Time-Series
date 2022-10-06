function [ts, wd] = genStaircase(steps, weights, fig)
% MakeFractal -- Create a deterministic discrete staircase (see Muzy, 1994)
%  Usage
%    Frac = MakeFractal(n,base,type,prob)
%  Inputs
%    steps    number of loops
%    weights  vector with weights
%               [.5, 0, .5] -> Monofractal Staircase
%               [.69, .46, -.46, .31] -> Multifractal Staircase
%  Outputs
%    ts       discrete time series (length ~ size(weights)^steps)
%    wd       weight distribution (length ~ size(weights)^steps)
%
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
       
        wd = wd.*weights(divs+1)';
        iTime = iTime+1;
    end
    
	ts = cumsum(wd);
    
    if fig == 1
        figure;
        plot(ts, 'k');
        title('Generalized Devil Staircase', 'FontSize', 24);
        xlabel('x', 'FontSize', 20);
        ylabel('f(x)', 'FontSize', 20);
    end
end

%  
%  Created by 
%       Gabriel Fornari
%  On 
%       18/05/2015 (dd/mm/yyyy)
%  
%  Based on
%       Part of Wavelab Version 850
%       Built Tue Jan  3 13:20:39 EST 2006
%       This is Copyrighted Material
%       For Copying permissions see COPYING.m
%       Comments? e-mail wavelab@stat.stanford.edu
