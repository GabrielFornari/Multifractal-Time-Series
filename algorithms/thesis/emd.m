function [IMFs] = emd(x, siftSteps, maxDif)
% Huang-Hilbert Transform (Empirical Mode Decomposition)
%   This transform decomposes a series x into a set of elemental 
%   signals called Intrinsic Mode Functions(IMFs). The IMFs are 
%   obtained by getting sifting the series with two envelopes: 
%   one upper and one lower. For more information see Huang et 
%   al(1999) and Huang et al(2003).
%
% Usage:
%   [IMFs] = emd(x, siftSteps, maxDif)
%
% Inputs
%   x           signal (time series)
%   siftSteps   optional; max number of steps(sift) to obtain an 
%               IMF
%   maxDif      optional; difference beetween global max and min of
%               the signal - used to stop sifting for more IMFs
%
% Outputs
%   IMFs - matrix N by M+1
%          N - length(X)
%          M - number of IMFs
%          (M+1) - residual (trend) of the series
%

    if nargin < 3
        maxDif = 0;
    end
    if nargin < 2
        siftSteps = 8;
    end
    
    % Read data and divide by std deviation
    xSize = length(x);
    xStd = std(x);
    x = x/xStd;

    % Evaluate TNM as total IMF number
    totalIMF = fix(log2(xSize))-1;
    
    % Find all IMFs
    nIMF = 1;
    while ~ismonotonic(x) && hasglobaldif(x, maxDif) ...
           && nIMF <= totalIMF
        hk = x;
        k = 1;
        while k <= siftSteps && ~isimf(hk)
            [maxValue, maxIndex] = findmax(hk);
            [minValue, minIndex] = findmin(hk);
            
            upp = spline(maxIndex, maxValue,1:xSize);
            low = spline(minIndex, minValue,1:xSize);
            
            hk = hk - (upp + low)/2;
            
            k = k+1;
        end      
        
        % Subtract IMF from data, use residual to find next IMF 
        x = x - hk;
        
        IMFs(:,nIMF) = hk(:);
        nIMF=nIMF+1;
    end

    % Put the residual (trend) in the last column  
    IMFs(:,end+1)=x(:);
    IMFs=IMFs*xStd;
end

function U = hasglobaldif(X, maxDiff)
% hasglobaldif   return if a series X has a global difference 
%                greater than max diff
%   A series X has its global difference calculated by the 
%   difference between its global maximum and global minimum.
%

    if max(X) - min(X) > maxDiff
        U = 1;
    else
        U = 0;
    end
end

function U = ismonotonic(X)
% ismonotonic   return if a series X is monotonic
%   A series is monotonic if it has only increasing or 
%   decreasing values, but not both. The series can also be 
%   constant.
%

    isAscending = all(diff(X) >= 0);
    isDescending = all(diff(X) <= 0);
    
    U = isAscending || isDescending;
end

function U = isimf(X)
% isimf    return if the series X is IMF
%   If the numbers of zero-crossings and extrema 
%   (local max + local min) are equal or differ by 1. This 
%   criterion was created by Huang et al.(1999, 2003).
%

    % Number os zero-crossings
    u1 = sum(X(1:end-1).*X(2:end) < 0);
    
    % Number os extrema (max + min)
    max = findmax(X);
    min = findmin(X);
    u2 = length(max)+length(min);
    
    if abs(u1-u2) > 1
        U = 0;
    else
        U = 1; 
    end
end

function [V, I] = findmax(X)
% findmax    local maximum of a series X
%   The peaks are obtained by getting the position which 
%   satisfies the following condition
%       x(i) >= x(i-1) AND x(i) >= x(i+1).
%
%   This function also changes the initial and end peaks to
%   improve the results of the funtion spline.
% 
%   [V, I] = findpeaks(X) returns the local maximum of a time 
%                         series
%       V - for VALUE of maximum
%       I - for INDEX of maximum
%

    xSize=length(X);
    
    max(1,1) = 1;
    max(1,2) = X(1);
    
    for iSize = 2:xSize-1
        if (X(iSize) >= X(iSize-1) && X(iSize) >= X(iSize+1))
            max(end+1,1) = iSize;
            max(end,2) = X(iSize);
        end
    end
    
    max(end+1,1) = xSize;
    max(end,2) = X(xSize);
    
    % Changes bounds to get a better spline
    if length(max(:, 1)) >= 4
        slope=(max(2,2)-max(3,2))/(max(2,1)-max(3,1));
        tmp=slope*(max(1,1)-max(2,1))+max(2,2);
        if tmp > max(1,2)
            max(1,2) = tmp;
        end
        slope=(max(end-1,2)-max(end-2,2))/...
              (max(end-1,1)-max(end-2,1));
        tmp=slope*(max(end,1)-max(end-1,1))+max(end-1,2);
        if tmp > max(end,2)
            max(end,2) = tmp;
        end
    end
    V = max(:, 2);
    I = max(:, 1);
end

function [V, I] = findmin(X)
% findmax    local minimum of a series X
%   The peaks are obtained by getting the position which 
%   satisfies the following condition
%       x(i) <= x(i-1) AND x(i) <= x(i+1).
%
%   This function also changes the initial and end peaks to
%   improve the results of the funtion spline.
% 
%   [V, I] = findpeaks(X) returns the local minimum of a time 
%                         series
%       V - for VALUE of minimum
%       I - for INDEX of minimum
%

    xSize=length(X);

    min(1,1) = 1;
    min(1,2) = X(1);
    
    for iSize = 2:xSize-1
        if X(iSize) <= X(iSize-1) && X(iSize) <= X(iSize+1)
            min(end+1,1) = iSize;
            min(end,2) = X(iSize);
        end
    end
    
    min(end+1, 1) = xSize;
    min(end, 2) = X(xSize);
    
    % Changes bounds to get a better spline
    if length(min(:, 1)) >= 4
        slope=(min(2,2)-min(3,2))/(min(2,1)-min(3,1));
        tmp=slope*(min(1,1)-min(2,1))+min(2,2);
        if tmp<min(1,2)
            min(1,2)=tmp;
        end
        slope=(min(end-1,2)-min(end-2,2))/...
              (min(end-1,1)-min(end-2,1));
        tmp=slope*(min(end,1)-min(end-1,1))+min(end-1,2);
        if tmp<min(end,2)
            min(end,2)=tmp;
        end
    end
    V = min(:, 2);
    I = min(:, 1);
end
