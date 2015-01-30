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
function [am, B, theta] = cirWaveCoefsFromSpread(M, s, k, varargin)
% Produces the curved (or short-crested) wave coefficient for curved
% incident waves from a cos squared spreading defined by s: 
%
%   B^2 = cos(1/2*theta).^(2*s)
%   
%   where B is the amplitude density spectrum.. 
% 
% Inputs:   M - coefficient limit (number of coefs is 2*M+1)
%           s - spreading factor 
%           k - the wave number
%           'RandPhase' - optional input, to create a random phase
%
% Ouputs:   am - circular wave coefficients
%           B - the directional spread
%           theta - theta points for B

opts = checkOptions({'RandPhase'}, varargin);
randPhase = opts(1);

% directional coefficients
Nthet = 4096;
dthet = 2*pi/Nthet;
theta = 0:dthet:(2*pi-dthet);
theta = -pi:dthet:(pi-dthet);

B2 = abs(cos(1/2*theta)).^(s);
A = trapz(theta,B2);
B = B2./A;

if (randPhase)
    B = B.*exp(2*pi*1i*rand(size(B)));
end

am = cirWaveCoefsFromFunc(B);
end

function [eps] = makePhase(theta)
    
    startEp = 2*pi*rand;
    
    N = length(theta);
    
    eps = zeros(1,N);
    eps(1) = startEp;
    ep = startEp;
    
    stepSi = 2*pi./N*100;
    
    for n = 2:N
        ep = ep + stepSi*randn;
        if (ep > 2*pi)
            ep = ep - 2*pi;
        elseif (ep < 0)
            ep = 2*pi - ep;
        end
        
        eps(n) = ep;
    end
    
end


