function [alpha, falpha, q, tauq] = mfwmm(x, wavelet, nvoice, oct, scale, fflag)				       									 
% Calculates the multifractal spectrum of a signal using Wavelet Transform
% Modulus Maxima method
% 
% Usage:
%   [alpha, falpha, q, tauq] = mfmm(x, wavelet, nvoice, oct, scale, figure)
%
% Inputs:
%   x        Signal (time series)
%   wavelet  String; can be: 'gaus', 'gaus1', 'gaus2', 'morl'
%   nvoice   number of nvoice/noctave, where
%   oct         noctave = floor(log2(n))-oct
%   scale       the minimum scale (scale *= 2)
%   figure   String; can be: 'wv', 'mm', 'tau', 'alp'
%                    or together: 'wvmmtau' (will plot 3 figures)
%
% Outputs:
%   alpha    Singularity Strength
%   falpha   Hausdorff Dimension
%   q        List of Exponents
%   tauq     Vector of Moments
%
    if nargin < 6
        fflag = 'none';
    end
    
    % Truncating the signal length to the biggest integer 2^n
    npts = 2^floor(log2(length(x)));
    w = x(1:npts);
    
    % Wavelet Transform
    [wta, scales] = cwt(w, wavelet, nvoice, oct, scale);
    
    %Wavelet Transform Modulus Maxima
    mma = wtmm(wta, 0.1);
    
    % Partition Function
    q = linspace(-20,20,41);
    za = thermoPart(wta, mma, q);
    
    % Power law             - Z(q, a) ~ a^tau(q)
    tauq = powerLaw(za, scales);
    
    % Singularity Spectrum  - D(h) = min(q*h - tau(q));
    [alpha, falpha] = legendreTransform(tauq, q);
    
    % Plots
    if ~isempty(strfind(fflag, 'wv'))  % Plot Wavelet Coeficients
        figure;
        hSurface = surf(linspace(0, 1, size(wta, 2)), log2(scales), wta);
        colormap(gray);
        set(hSurface, 'EdgeAlpha', 0);
        title('Wavelet Transform', 'FontSize', 16);
        xlabel('x', 'FontSize', 14);
        ylabel('log_{2}(a)', 'FontSize', 14);
    end
    if ~isempty(strfind(fflag, 'mm'))  % Plot Wavelet Modulus Maxima
        figure;
        hSurface = surf(linspace(0, 1, size(wta, 2)), log2(scales), mma);
        colormap(gray);
        set(hSurface, 'EdgeAlpha', 0);
        title('Modulus Maxima Map', 'FontSize', 16);
        xlabel('x', 'FontSize', 14);
        ylabel('log_{2}(a)', 'FontSize', 14);
    end
    if ~isempty(strfind(fflag, 'tau'))  % Plot tau vs. q
        figure;
        plot(q, tauq, 'ko');
        title('\tau(q) versus q', 'FontSize', 16);
        xlabel('q', 'FontSize', 14);
        ylabel('\tau(q)', 'FontSize', 14);
    end
    if ~isempty(strfind(fflag, 'alp'))  % Plot falpha vs. alpha
        figure;
        plot(alpha, falpha, 'kx');
        title('Singularity Spectrum - D(h) versus h', 'FontSize', 16);
        xlabel('h', 'FontSize', 14);
        ylabel('D(h)', 'FontSize', 14);
    end
end
%**********************************************************
%  Continuous Wavelet Transform
%  
%  Usage
%    wt = cwt(x, wavelet, nvoice, oct, scale)
%  Inputs
%    x        signal, dyadic length n=2^J, real-valued
%    wavelet  string 'gaus', 'gaus1', 'gaus2', 'morl'
%    nvoice   number of voices/octave
%    oct      Default=2
%    scale    Default=4
%  Outputs
%    wt       Wavelet Coeficients (matrix(nscale, time))
%         where,
%               time = length(x)
%               nscale = nvoice.*noctave (from low to high frequencies)
%
function [cwt, scales] = cwt(x, wavelet, nvoice, oct, scale)
    if nargin<5
		oct = 2;
		scale = 4;
    end
    
	n = length(x);
	xhat = fft(x);
	xi = [(0:(n/2)) (((-n/2)+1):-1)]' .* (2*pi/n);
    
	noctave = floor(log2(n))-oct;
	nscale  = nvoice*noctave;
    scales = zeros(nscale, 1);
    
    cwt = zeros(nscale, n);
    omega0 = 5;
	kscale  = 1;

	for jo = 1:noctave
        for jv = 1:nvoice
            qscale = scale*(2^(jv/nvoice));
            omega =  n .* xi ./ qscale;
            if strcmp(wavelet,'gaus')
                window = exp(-omega.^2 ./2);
            elseif strcmp(wavelet,'gaus1')
                window = 1i.*omega.*exp(-omega.^2 ./2);
            elseif strcmp(wavelet,'gaus2')
                window = (omega.^2) .* exp(-omega.^2 ./2);
            elseif strcmp(wavelet,'morl')
                window = exp(-(omega - omega0).^2 ./2) - exp(-(omega.^2 + omega0.^2)/2);
            end
            
            scales(kscale) = qscale;
            what = window.*xhat;
            w    = ifft(what);
            cwt(nscale+1-kscale, 1:n) = real(w);
            kscale = kscale+1;
        end
		scale = scale.*2;
    end
end
%
%  Modified by 
%       Gabriel Fornari
%  On 
%       09/03/2015 (dd/mm/yyyy)
%  
%  Based on
%       Originally created for WaveLab.701.
%
%       Modified by Maureen Clerc and Jerome Kalifa, 1997
%       clerc@cmapx.polytechnique.fr, kalifa@cmapx.polytechnique.fr   
%   
%       Part of WaveLab Version 802
%       Built Sunday, October 3, 1999 8:52:27 AM
%       This is Copyrighted Material
%       For Copying permissions see COPYING.m
%       Comments? e-mail wavelab@stat.stanford.edu
%**********************************************************

%**********************************************************
%  Modulus Maxima of wavelet transform
%
%  Usage
%    mmap = wtmm(wt, threshold)
%  Inputs
%    wt         Wavelet Transform     (matrix(scales, time))
%    threshold  Minimum Peak Height in percent (default = 0.5)
%  Outputs
%    mmap       Binary indicating max (matrix - n, m | n = scales; m = time)
%
function mmap = wtmm(wt, threshold)
    if nargin == 1
        threshold = .5;
    end

    nScale = size(wt, 1);
    nTime = size(wt, 2);
	mm = zeros(size(wt));
    mmap = zeros(size(wt));
    
    wt = abs(wt);
    
    for iScale=1:nScale
        for iTime=2:nTime-1;
            if ((wt(iScale, iTime)>=wt(iScale, iTime-1)) && (wt(iScale, iTime)>wt(iScale, iTime+1)))
                mm(iScale, iTime)=wt(iScale, iTime);
                mmap(iScale, iTime)=1;
            end
        end
    end
    % Remove threshold values
    for iScale=1:nScale
        temp = max(mm(iScale, :));
        mm(iScale, :)=mm(iScale, :)/temp;
        for iTime=2:nTime-1;
            if mm(iScale, iTime)<threshold
                mmap(iScale, iTime)=0;
            end
        end
    end
end
%  
%  Created by 
%       Gabriel Fornari
%  On 
%       03/03/2015 (dd/mm/yyyy)
%**********************************************************    

%**********************************************************
%  Thermodynamic Partition Function
%
%  Usage
%    z = thermoPart(cw, mm, q)
%  Inputs
%    cw    Wavelet Coeficients  (matrix(scales, time)) --> See "cwt"
%    mm    Binary Map           (matrix(scales, time)) --> See "mmwt"
%    q     List of Exponents    (default = linspace(-20,20,61))
%  Outputs
%    z     matrix nexp by nscale of z(a, q)
%
%  Description
%    z(q, a) = sum(|CWT(a,b(i))|^q)
%           where b = (b(i)) is a list of local maxima  
%
function z = thermoPart(cw,mm,q)
	if nargin < 3
	   q = linspace(-20,20,61);
    end
	
    nScale = size(cw, 1);
    z = zeros(nScale, length(q));
    
	for iScale=1:nScale
	    j = find(mm(iScale, :));
		if ~isempty(j)
            c = abs(cw(iScale, j));
            for i=1:length(q)
                z(iScale, i) = sum(c.^q(i));
            end
		else
			z(iScale, :) = eps.^q(:);
        end
    end
end
%  
%  Modified by 
%       Gabriel Fornari (gabriel.fornari@inpe.br)  
%  At 
%       09/03/2015 (dd/mm/yyyy)
%  
%  Based on
%       Wavelab Version 850
%       Built Tue Jan  3 13:20:39 EST 2006
%       This is Copyrighted Material
%       For Copying permissions see COPYING.m
%       Comments? e-mail wavelab@stat.stanford.edu
%**********************************************************

%**********************************************************
%  Power Law behavior -> Z(a, q) ~ a^(tau(q))
%  
%  Usage
%    tau = powerLaw(z, scales)
%  Inputs
%    z         Values of partition function (matrix (q, scales))
%    scales    List of Scales
%  Outputs
%    tau       Vector of moments (length(tau) = length(scale))
%
function tau = powerLaw(z, scales)
    nq = size(z, 2);
    tau = zeros(nq, 1);
    
    for iq = 1:nq
        coef = polyfit(log(scales), log(z(:, iq)), 1);
        tau(iq) = coef(1);
    end
end
%  
%  Modified by 
%       Gabriel Fornari (gabriel.fornari@inpe.br)  
%  On 
%       09/03/2015 (dd/mm/yyyy)
%  
%  Based on
%       Wavelab Version 850
%       Built Tue Jan  3 13:20:39 EST 2006
%       This is Copyrighted Material
%       For Copying permissions see COPYING.m
%       Comments? e-mail wavelab@stat.stanford.edu
%**********************************************************

%**********************************************************
%  Singularity Spectrum of Local Scaling Exponents
%  
%  Usage
%    [alpha, falpha] = singSpec(tau, q)
%  Inputs
%    tau       vector (length(tau) = length(q))
%    q         list of exponents
%  Outputs
%    alpha     the same of input but without falpha negative values
%    falpha    vector without negative fractal dimensions
%
%  Description
%    falpha(alpha) = min[for each q](alpha*q - tau(q))
%
function [alpha, falpha] = legendreTransform(tau, q)
    alpha = diff(tau);
    for iAlpha = 1:length(tau)-1
        alpha(iAlpha) = alpha(iAlpha)/(q(iAlpha+1)-q(iAlpha));
    end
    falpha = (q(1:end-1).*alpha') - tau(1:end-1)';
    
    % Remove negative values
    alpha(sign(falpha)==-1) = [];
    falpha(sign(falpha)==-1) = [];
end
%  
%  Created by 
%       Gabriel Fornari (gabriel.fornari@inpe.br)  
%  On 
%       02/04/2015 (dd/mm/yyyy)
%**********************************************************
