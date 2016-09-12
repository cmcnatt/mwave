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
function [f, a, ep, off, R, S, freqs] = cosFit(x, y, varargin)

L = length(y);

off0 = mean(y);

Y = fft(y-off0)/L;
S = 2*Y(1:floor(L/2)+1);
freqs = (0:floor(L/2))/(x(2)-x(1))/L;

[a0, ind] = max(abs(S));
f0 = freqs(ind);
ep0 = angle(S(ind));
%ep0 = wrapTo2Pi(ep0 + pi);

% Refine the fit (this does funny things...)
% sol = fminsearch(@(param) sum((y(:) - param(4) - param(2)*cos(2*pi*param(1)*x(:) + param(3))).^2), [f0 a0 ep0 off0]);

% try it on just 50 cycles.
% T0 = 1/f0;
% [~, iStop] = min(abs(x - 50*T0));
% 
% 
% sol = fminsearch(@(param) sum((y(1:iStop) - param(4) - param(2)*cos(2*pi*param(1)*x(1:iStop) + param(3))).^2), [f0 a0 ep0 off0]);
% 
% f = sol(1);
% a = sol(2);
% ep = sol(3);
% off = sol(4);

% use the curve fitting toolbox;
% if (isrow(x))
%     x = x.';
% end
% if (isrow(y))
%     y = y.';
% end
% type = fittype('a*cos(2*pi*f*x + ep) + off');
% cosFit = fit(x, y, type, 'StartPoint', [a0, ep0, f0, off0]);
% %cosFit = fit(x, y, type, 'StartPoint', [a0, ep0, f0, off0], 'Lower', [f0-0.1*f0, a0-0.5*a0, 0, -Inf], 'Upper', [f0+0.1*f0, a0 + 0.5*a0, 2*pi, Inf]);
% 
% f = cosFit.f;
% a = cosFit.a;
% ep = cosFit.ep;
% off = cosFit.off;
% 
% f = f0;
% a = a0;
% ep = ep0;
% off = off0;
% 
% if (a < 0)
%     a = -a;
%     ep = ep + pi;
%     if (ep >= 2*pi)
%         ep = ep - 2*pi;
%     end
% end

% sol = fminsearch(@(param) sum((y(:) - off0 - a0*cos(2*pi*param(1)*x(:) + ep0)).^2), f0);
% 
% f = sol(1);
% a = a0;
% ep = ep0;
% off = off0;


f = f0;
a = a0;
ep = ep0;
off = off0;

y0 = a*cos(2*pi*f*x + ep) + off;

R = sum((y - y0).^2);

end

