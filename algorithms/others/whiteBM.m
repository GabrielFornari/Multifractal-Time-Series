function ts = whiteBM(n, fig)
% Brownian Motion
%   Header missing
    
    mean = 0;
    variance = 1;
    
    ts = zeros(n, 1);
    
    for i=2:n
        ts(i) = ts(i-1)+normrnd(mean, sqrt(variance))/(n-1);
    end
    
    if fig == 1
        figure;
        plot(ts, 'k');
        title('Brownian Motion', 'FontSize', 16);
        xlabel('x', 'FontSize', 14);
        ylabel('f(x)', 'FontSize', 14);
    end
end
