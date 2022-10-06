function x = gpmodel(steps, p, l, fflag, n)
% Creates a multifractal series based on the p-model.
%
% Usage
%   [x] = pmodelSignal(steps, p, l, fflag, n);
%
% Example
%   [x] = pmodelSignal(11, .3, .5);
%
% Inputs
%   steps    number of loops
%   p        probability (beetween 0 and 1)
%   l        segment length (beetween 0 and 1)
%   fflag    1|0 - for output plot
%   n        sample length (optional)
%
% Outputs
%   x       multifractal signal (length = 2^steps OR length = n)
%
    
    if nargin < 5
        fflag = 0;
    end
    if nargin < 6
        n = 2^steps;
    end
    
    x = ones(n, 1);
    p1 = p;
    p2 = 1-p;
    
    l1 = l;
    l2 = 1-l1;
    
    for iStep = 0:steps-1
        nSgm = 2^iStep;
        sgm = round(n/nSgm);
        
        for iSgm = 1:nSgm
            switch randi([1 4]);
                case 1
                    x(((iSgm-1)*sgm)+1:((iSgm-1)*sgm)+sgm*l1) = x(((iSgm-1)*sgm)+1:((iSgm-1)*sgm)+sgm*l1)*2*p1;
                    x(((iSgm-1)*sgm)+1+sgm*l2:(iSgm*sgm)) = x(((iSgm-1)*sgm)+1+sgm*l2:(iSgm*sgm))*2*p2;
                case 2
                    x(((iSgm-1)*sgm)+1:((iSgm-1)*sgm)+sgm*l1) = x(((iSgm-1)*sgm)+1:((iSgm-1)*sgm)+sgm*l1)*2*p2;
                    x(((iSgm-1)*sgm)+1+sgm*l2:(iSgm*sgm)) = x(((iSgm-1)*sgm)+1+sgm*l2:(iSgm*sgm))*2*p1;
                case 3
                    x(((iSgm-1)*sgm)+1:((iSgm-1)*sgm)+sgm*l2) = x(((iSgm-1)*sgm)+1:((iSgm-1)*sgm)+sgm*l2)*2*p1;
                    x(((iSgm-1)*sgm)+1+sgm*l1:(iSgm*sgm)) = x(((iSgm-1)*sgm)+1+sgm*l1:(iSgm*sgm))*2*p2;
                case 4
                    x(((iSgm-1)*sgm)+1:((iSgm-1)*sgm)+sgm*l2) = x(((iSgm-1)*sgm)+1:((iSgm-1)*sgm)+sgm*l2)*2*p2;
                    x(((iSgm-1)*sgm)+1+sgm*l1:(iSgm*sgm)) = x(((iSgm-1)*sgm)+1+sgm*l1:(iSgm*sgm))*2*p1;
            end
        end
    end
    
    if fflag
        figure;
        plot(x, 'k');
    end
end
