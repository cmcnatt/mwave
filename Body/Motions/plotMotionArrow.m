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
function [] = plotMotionArrow(type, cent, len, varargin)

linwid = 1;
arrlen = 0.1*len*cosd(45);

isrot = true;

if (strcmp(type,'surge') || strcmp(type,'sway') || strcmp(type,'heave'))
    isrot = false;
end

if (isrot)
    r = varargin{1};
    
    thet1 = cent - 0.5*len;
    thet2 = cent + 0.5*len;

    theta = thet1:(thet2-thet1)/100:thet2;
end

switch type
    case 'surge'
        norm = [1 0 0];

        a2t = [-arrlen, 0, arrlen];
        a2b = [-arrlen, 0, -arrlen]; 

        a1t = [arrlen, 0, arrlen];
        a1b = [arrlen, 0, -arrlen];
    case 'sway'
        norm = [0 1 0];

        a2t = [0, -arrlen, -arrlen];
        a2b = [0, -arrlen, arrlen]; 

        a1t = [0, arrlen, arrlen];
        a1b = [0, arrlen, -arrlen];
    case 'heave'
        norm = [0 0 1];

        a2t = [-arrlen, 0, -arrlen];
        a2b = [arrlen, 0, -arrlen]; 

        a1t = [arrlen, 0, arrlen];
        a1b = [-arrlen, 0, arrlen];
    case 'roll'
        xr = zeros(size(theta));
        yr = r*cos(theta);
        zr = r*sin(theta);
    case 'pitch'
        
        arrlen = 2*arrlen;
        xr = r*cos(theta);
        yr = zeros(size(xr));
        zr = r*sin(theta);
        
        p1 = r*[cos(thet1) 0 sin(thet1)];
        p2 = r*[cos(thet2) 0 sin(thet2)];
        
        ang2 = -(thet2-pi/2-pi);
        
        R2 = [cos(ang2) 0 sin(ang2); 0 1 0; -sin(ang2) 0 cos(ang2)];
        
        a2t = (R2*[-arrlen, 0, arrlen]')';
        a2b = (R2*[-arrlen, 0, -arrlen]')'; 

        ang1 = -(thet1+pi/2);
        R1 = [cos(ang1) 0 sin(ang1); 0 1 0; -sin(ang1) 0 cos(ang1)];
        a1t = (R1*[arrlen, 0, arrlen]')';
        a1b = (R1*[arrlen, 0, -arrlen]')';
        
    case 'yaw'
        xr = r*cos(theta);
        yr = r*sin(theta);
        zr = zeros(size(xr));
end

if (~isrot)
    p1 = cent - 0.5*len*norm;
    p2 = cent + 0.5*len*norm;

    plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'k', 'linewidth', linwid);    
else
    plot3(xr, yr, zr, 'k', 'linewidth', linwid);
end

p2t = p2 + a2t;
p2b = p2 + a2b;

p1t = p1 + a1t;
p1b = p1 + a1b;

hold on;
plot3([p2t(1) p2(1) p2b(1)], [p2t(2) p2(2) p2b(2)], [p2t(3) p2(3) p2b(3)], 'k', 'linewidth', linwid);
plot3([p1t(1) p1(1) p1b(1)], [p1t(2) p1(2) p1b(2)], [p1t(3) p1(3) p1b(3)], 'k', 'linewidth', linwid);


