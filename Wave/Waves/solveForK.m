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
function [k] = solveForK(omega, h, g, L)
% k = SolveForK(omega, h, g, L)
% solves for the wavenumber as a function of radial wave frequency (omega)
% and depth (h).
%
% g is optional. [] will use default value.
% 
% Optional argument L which request the evanescent wave numbers as well 
% where L is an integer which is the number of requested evanescent wave 
% numbers.

if nargin < 3
    g = [];
end

if isempty(g)
    g = 9.80665;
end

if nargin < 4
    L = [];
end

evan = false;
if ~isempty(L)
    evan = true;
end

if (isempty(omega) || isempty(h) || isempty(g))
    error('Empty parameters in solveForK.');
end

if (h <= 0)
    error('The depth must be greater than 0');
elseif (h == Inf)
    k = omega.^2/g;
else
    tol = 1e-10;
    k = omega.^2/g;
    const = h.*k;

    k0h = k.*h;
    tanhk0h = tanh(k0h);
    f0 = const - k0h.*tanhk0h;

    i = 0;

    while (any(f0 > tol))
        i = i+1;
        k0h = k.*h;
        tanhk0h = tanh(k0h);
        f0 = const - k0h.*tanhk0h;
        m = k0h.*tanhk0h.^2 - tanhk0h - k0h;
        kh = k0h - f0./m;
        k = kh./h;
    end
end

if (evan)
    if (~isInt(L))
        error('If evanescent modes are requested, the argument after ''Evanescent'' must be an integer, which is the reqeusted number of evanescent modes');
    end

    kL = zeros(1,L);

    if (~isinf(h))
        % Bisection method - always converges
        for n = 1:L
            tol = 1e-7;
            const = omega.^2.*h/g;

            ul = n*pi;
            ll = n*pi - pi/2;

            f0 = 1;
            i = 0;

            while (abs(f0) > tol)
                i = i+1;
                kh = (ul + ll)/2;
                tankh = tan(kh);
                f0 = const + kh.*tankh;

                if (f0 < 0)
                    ll = kh;
                else
                    ul = kh;
                end

                if (i > 1e3)
                    break;
                end
            end

            kL(n) = kh./h;
        end

        %{
        % Newton-Raphson method
        for n = 1:N
            tol = 1e-8;

            const = omega.^2.*h/g;

            k0h = pi*(n - 1/2 + 0.01/h);
            kn = k0h./h;
            tank0h = tan(k0h);
            f0 = const + k0h.*tank0h;

            i = 0;

            while (any(abs(f0) > tol))
                i = i+1;
                k0h = kn.*h;
                tank0h = tan(k0h);
                f0 = const + k0h.*tank0h;
                m = tank0h + k0h + k0h.*tank0h.^2;
                kh = k0h - f0./m;
                kn = kh./h;
            end

            kN(n) = kn;
        end
        %}
    end

    k0 = k;
    k = zeros(1, L+1);
    k(1) = k0;
    k(2:end) = kL;
end
