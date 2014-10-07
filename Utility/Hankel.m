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
function [H] = Hankel(m, z)

if(~isInt(m))
    error('m must be an integer');
end

if (m == 0 || m == 1)
    H = besselh(m, 2, z);
else
    H0 = besselh(0, 2, z);
    H1 = besselh(1, 2, z);
    
    for n = 2:m
        H = (2*(n-1)./z).*H1 - H0;
        H0 = H1;
        H1 = H;
    end
end
