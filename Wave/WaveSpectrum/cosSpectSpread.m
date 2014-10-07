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
function [g] = cosSpectSpread(s, thetac, theta)
% Reference: Holthuijsen - Waves is Oceanic Waters - eq 6.3.25
% s - spreading width parameter
% thetac - center direction (rad)
% theta - input directions (rad)
%  The intregral of the spreading function is 1.

if (s <= 0)
    error('pow must be greater than 0');
elseif (mod(s, 2) ~= 0)
    error('pow must be an even number');
end

if (abs(thetac) > pi)
    error('thetac must be between -pi and pi');
end

A2 = gamma(s+1)/(gamma(s+1/2)*2*sqrt(pi));
g = A2*cos(1/2*(theta - thetac)).^(2*s);

delTh = computeDelta(theta);
gArea = sum(g.*delTh);

g = g./gArea;


