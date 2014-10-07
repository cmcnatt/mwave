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
function [inter] = segmentsIntersect(pnt11, pnt12, pnt21, pnt22)

if (length(pnt11) ~= 2)
    error('each point must contain 2 values: {x,y}');
end

if (length(pnt12) ~= 2)
    error('each point must contain 2 values: {x,y}');
end

if (length(pnt21) ~= 2)
    error('each point must contain 2 values: {x,y}');
end

if (length(pnt22) ~= 2)
    error('each point must contain 2 values: {x,y}');
end

x11 = pnt11(1);
x12 = pnt12(1);
x21 = pnt21(1);
x22 = pnt22(1);

m1 = (pnt12(2) - pnt11(2))/(x12 - x11);
m2 = (pnt22(2) - pnt21(2))/(x22 - x21);

b1 = pnt11(2) - m1*x11;
b2 = pnt21(2) - m2*x21;

if (m1 == m2)
    if (b1 == b2)
        inter = true;
    else
        inter = false;
    end
else
    b1 = pnt11(2) - m1*x11;
    b2 = pnt21(2) - m2*x21;
    
    xint = (b2 - b1)/(m1 - m2);
    
    if (x12 >= x11)
        if ((xint <= x11) || (xint > x12))
            inter = false;
            return
        end
    else
        if ((xint < x12) || (xint >= x11))
            inter = false;
            return
        end
    end
    
    if (x22 >= x21)
        if ((xint <= x21) || (xint > x22))
            inter = false;
            return
        end
    else
        if ((xint < x22) || (xint >= x21))
            inter = false;
            return
        end
    end
    
    inter = true;
end