function ts = fracBM(n, h, fig)
% Brownian Motion
%   Header missing
    
    mean = 0;
    variance = 1;
    beta = 2*h+1
    
    ts = zeros(n, 1);
    
    x = zeros(n/2, 1);
    
    for i=1:n/2
        rad = i^(-beta/2)*normrnd(mean, sqrt(variance))/(n-1);
        phase = 2*pi*rand;
        
        x(i) = rad*(cos(phase)+1j*sin(phase));
    end
    
    ts = real(ifft(x));
    
    if fig == 1
        figure;
        plot(ts, 'k');
        title('Brownian Motion', 'FontSize', 16);
        xlabel('x', 'FontSize', 14);
        ylabel('f(x)', 'FontSize', 14);
    end
end
