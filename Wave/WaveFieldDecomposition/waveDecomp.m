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
function [A, eta0] = waveDecomp(k, r0, theta, z, h, eta, L, M, varargin)
%   Produces the complex coefficients needed for a circular wave field from
%   a non-dimensional pressure, p/(rho*g), which is the input eta0.  eta0 
%   is the value at a distance, r, as a function of direction, theta, and 
%   depth, z, for a wave field with wave number, k.
%
%   The number of depths determines how many evanescent wave mode
%   coefficients are computed.  If only a single depth is present, only the
%   propagating modes are computed and the number of evanescent modes, M, 
%   is length(z) - 1.  
%
%   In the cicular direction, the number of coefficients is specified by
%   the input parameter, M.  In the circular direction, the summation goes 
%   from -M to M and so the number of the coefficients is 2M+1.  But if the
%   optional argument 'TrimCoeffs' is provided the bound of the summation 
%   is trimmed to some lower value of M' if and only if the magnitude of 
%   the coefficients beyond M' are infintesimally small.  That is, for some 
%   cases only M' coefficients are required even though the user specified
%   M.
%
%   The output matrix of coefficients is of size (L + 1) x (2M + 1). The 
%   first row contains the coefficients of the propagating modes and each 
%   subsequent row contains the coefficients of the subsequent evanescent 
%   mode. The columns are the coefficients for the circular modes for the
%   propagating mode then each evanescent mode.
%
%   Optional arguments:
%     - 'SigFigCutoff', val - val is the number of significant figures to 
%       keep. Cut off all values that insignificant.
%       [A, eta0] = waveDecomp(k, r0, theta, z, h, eta, L, M, ...
%               'SigFigCutoff', val)
%     - 'Incident' - computes the coefficients for an incident cylindircal
%       wave rather than an outwardly radiating one

[opts, args] = checkOptions({{'SigFigCutoff', 1}, {'Incident'}, {'Round', 1}}, varargin);
useSigFig = opts(1);
if (useSigFig)
    sigFig = args{1};
else
    sigFig = -1;
end

isInc = opts(2);
if (isInc)
    L = 0;
end

round2 = 0;
rnd = false;
if (opts(3))
    rnd = true;
    round2 = args{3};
end

[Nz, Ntheta] = size(eta);

if (length(theta) ~= Ntheta)
    error('The number of eta colums must be the same as the number of theta points.');
end

if (length(z) ~= Nz)
    error('The number of eta rows must be the same as the number of z points.');
end

if (z(1) ~= 0)
    error('The first z point must be zero');
end

for n = 1:Ntheta
    if (theta(n) < 0)
        error('theta values must be between 0 and 2*pi inclusive');
    end
end

% integrate over the depth first

k0 = k(1);
kl = k(2:end);

if (length(z) == 1)
    L = 0;
    etaL = eta;
else
    if ((z(end) + h) > 1e-4)
        error('The last z point must be at the bottom, z = -h');
    end

    z = flipud(z);
    etaz = flipud(eta);

    etaL = zeros(L+1,Ntheta);

    % l = 0
    coshkz = cosh(k0*(z + h));
    const0 = 2*cosh(k0*h)/(h*(1+sinh(2*k0*h)/(2*k0*h)));

    for n = 1:Ntheta
        etaL(1, n) = const0*trapz(z, etaz(:,n).*coshkz);
    end

    constl = 2./(h*(1 + sin(2*kl*h)./(2*kl*h)));
    cosklz = cos((z + h)*kl);

    for l = 1:L
        for n = 1:Ntheta
            etaL(l+1, n) = constl(l)*trapz(z, etaz(:,n).*cosklz(:,l));
        end
    end
end

% then find the circular coefficients
% l = 0 - progressive coefficients
if (~isInc)
    H0 = besselh(0, 2, k0*r0);
    Hm = besselh(1:M, 2, k0*r0);
    Hnm = (-1).^(1:M).*Hm;
else
    H0 = besselj(0, k0*r0);
    Hm = besselj(1:M, k0*r0);
    Hnm = (-1).^(1:M).*Hm;
end

X = 1/Ntheta*fft(etaL(1,:));

if(rnd)
    X = round(X./round2).*round2;
end

A0 = zeros(1, 2*M+1);
A0(M+1) =  X(1)/H0;
A0(M+2:2*M+1) = X(2:M+1)./Hm;
A0(M:-1:1) = X(Ntheta:-1:Ntheta-M+1)./Hnm;

if(rnd)
    A0 = round(A0./round2).*round2;
end

Mold = M;

if (sigFig > 0)
    maxVal = max(abs(etaL(1,:)));
    expV = floor(log10(maxVal));
    cutOffVal = 10^(expV-sigFig);
    
    trimM = M;
    
    for m = M:-1:1
        Am = abs(A0(M+1+m));
        Anm = abs(A0(M+1-m));

        if (all(cutOffVal >= Am) && all(cutOffVal >= Anm))
            trimM = m-1;
        else
            break;
        end
    end
    
    M = trimM;
end

% the output coeffcient matrix
A = zeros(L + 1, 2*M + 1);

A(1,:) = A0(Mold+1-M:Mold+1+M);

% the evanescent modes
for l = 1:L
    X = 1/Ntheta*fft(etaL(l+1,:));
    
    K0 = besselk(0, kl(l)*r0);
    Km = besselk(1:M, kl(l)*r0);
    
    A(l+1,M+1) = X(1)/K0;
    A(l+1,M+2:2*M+1) = X(2:M+1)./Km;
    A(l+1,M:-1:1) = X(Ntheta:-1:Ntheta-M+1)./Km;
end

eta0 = etaL(1,:);

end