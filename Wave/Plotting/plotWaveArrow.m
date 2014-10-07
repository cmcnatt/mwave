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
function [] = plotWaveArrow(start, beta, len)

linwid = 1;

p1 = start;

p2 = p1 + [len*cos(beta), len*sin(beta), 0];

arrlen = 0.1*len*cosd(45);
p2t = p2 + [-arrlen, 0, arrlen];
p2b = p2 + [-arrlen, 0, -arrlen]; 

plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'linewidth', linwid);
hold on;

plot3([p2t(1) p2(1) p2b(1)], [p2t(2) p2(2) p2b(2)], [p2t(3) p2(3) p2b(3)], 'linewidth', linwid);

wlen = 0.75*len;

val = [0.3; 0.5; 0.7];

wxpts = ones(3,1)*p1 + val*[len*cos(beta), len*sin(beta), 0];

wxpts1 = wxpts + 0.5*wlen*ones(3,1)*[-sin(beta), cos(beta), 0];
wxpts2 = wxpts + 0.5*wlen*ones(3,1)*[sin(beta), -cos(beta), 0];

plot3([wxpts1(1,1), wxpts2(1,1)], [wxpts1(1,2), wxpts2(1,2)], [wxpts1(1,3), wxpts2(1,3)],'linewidth', linwid)
plot3([wxpts1(2,1), wxpts2(2,1)], [wxpts1(2,2), wxpts2(2,2)], [wxpts1(2,3), wxpts2(2,3)],'linewidth', linwid)
plot3([wxpts1(3,1), wxpts2(3,1)], [wxpts1(3,2), wxpts2(3,2)], [wxpts1(3,3), wxpts2(3,3)],'linewidth', linwid)

