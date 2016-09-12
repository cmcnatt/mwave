%{ 
mwave - A water wave and wave energy converter computation package 
Copyright (C) 2014  Cameron McNatt

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Contributors:
    C. McNatt
%}
function [f, a, ep, off, S, freqs, a2, ep2] = regularWaveAnalysis(sampleFreq, sig, varargin)

if iscolumn(sig)
    sig = sig.';
end
N = length(sig);
dt = 1/sampleFreq;
t = 0:dt:(N-1)*dt;


[opts, args] = checkOptions({{'StartTime', 1}, {'StopTime', 1}, ...
    {'ExpectFreq', 1} {'Window'}, {'Range', 1}, {'AdjPhase'}}, varargin);

startT = 0;
startI = 1;
if (opts(1))
    startT = args{1};
    [~, startI] = min(abs(t - startT));
    startT = t(startI);
end

stopT = t(end);
stopI = N;
if (opts(2))
    stopT = args{2};
    [~, stopI] = min(abs(t - stopT));
    stopT = t(stopI);
end

sig2 = sig(startI:stopI);
t2 = t(startI:stopI);

useWind = opts(4);
range = [];
if opts(5)
    range = args{5};
end
adjPhase = opts(6);


% guess 1
if (opts(3))
    f1 = args{3};
    off = mean(sig2);
else
    [f1, a1, ep1, off, ~, S1, freqs1] = cosFit(t2, sig2);
end

% use guess 1 to get an exact number of cycles in the window.
T = 1/f1;

DeltaT = stopT - startT;
Ncyl = floor(DeltaT/T);

% guess 2
if (Ncyl > 1)
    stopT2 = startT + T*Ncyl;
    [~, stopI2] = min(abs(t - stopT2));
    
    % new bit
    sigSval = sig(startI);
    isUp = false;
    if (sig(startI+1) > sigSval)
        isUp = true;
    end
    
    cycIndLen = floor(T*sampleFreq);
    
    sigEnd = sig(stopI2-cycIndLen:stopI2);
    minVal = 1000;
    minValInd = 0;
    for n = 2:length(sigEnd)
        val = abs(sigEnd(n) - sigSval);
        
        if (sigEnd(n) > sigEnd(n-1))
            if (isUp)
                if (val < minVal)
                    minVal = val;
                    minValInd = n;
                end
            end
        else
            if (~isUp)
                if (val < minVal)
                    minVal = val;
                    minValInd = n;
                end
            end
        end
    end
    
    stopI2 = stopI2 - length(sigEnd) + minValInd - 1;
    % end new bit
    
    sig3 = sig(startI:stopI2) - off;
    t3 = t(startI:stopI2);
    if (useWind)
        M = length(sig3); 
        w = tukeywin(M, 0.1);
        %w = hanning(M); 
        sig3 = w'.*sig3;
        figure;
        plot(t3, sig3);
    end
    [f, a, ep, ~, ~, S, freqs] = cosFit(t3, sig3);
else
    f = f1;
    a = a1;
    ep = ep1;
    S = S1;
    freqs = freqs1;
end

a2 = 0;
ep2 = 0;
[maxS, ind] = max(abs(S));
ind2 = 2*(ind - 1) + 1;

if ~isempty(range)
    istart = indexOf(freqs, freqs(ind)-range/2);
    istop = indexOf(freqs, freqs(ind)+range/2);
    
    Snew = sum(S(istart:istop));
    
    istart = indexOf(freqs, freqs(ind2)-range/2);
    istop = indexOf(freqs, freqs(ind2)+range/2);
    Snew2 = sum(S(istart:istop));
    
    a = abs(Snew);
    ep = angle(Snew);
    
    a2 = abs(Snew2);
    ep2 = angle(Snew2);
else
    % second order coefficients
    if (ind2 < length(S))
        a2 = abs(S(ind2));
        ep2 = angle(S(ind2));
    end
end

if (adjPhase)
    x = t3;
    y = sig3;
    delEps = linspace(0,2*pi,101);
    
    minR = 1e10;
    iminR = 1;
    
    for n = 1:length(delEps)
        delEp = delEps(n);
        y0 = a*cos(2*pi*f*x + ep + delEp) + a2*cos(2*2*pi*f*x + ep2 + delEp);
        %y0 = a*cos(2*pi*f*x + ep + delEp) + off;
        R = sum((y - y0).^2);
        if (R < minR)
            minR = R;
            iminR = n;
        end
    end
    ep = ep + delEps(iminR);
    ep2 = ep2 + delEps(iminR);
end

