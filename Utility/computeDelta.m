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
function [dx] = computeDelta(x)

[row col] = size(x);

if (row > 1 && col > 1)
    error('Input must be a vector');
end

n = length(x);

midx = (x(2:n) + x(1:n-1))./2;
midx1 = 2*x(1) - midx(1);
midxn = 2*x(n) - midx(n-1);

if (col > 1)
    midx = [midx1, midx, midxn];
else
    midx = [midx1; midx; midxn];
end

dx = midx(2:n+1) - midx(1:n);

end

