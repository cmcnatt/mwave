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
function [xo, yo] = resampleEven(xi, yi, N)

dx = diff(xi);
dy = diff(yi);

len = sqrt(dx.^2 + dy.^2);
del = cumsum([0 len]);

delta = sum(len)/(N-1);

xii = 0:delta:(N-1)*delta;

xo = interp1(del, xi, xii, 'linear', 'extrap');
yo = interp1(del, yi, xii, 'linear', 'extrap');