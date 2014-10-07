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
function [geo] = makePanel_spheroidEndHorCyl(length, radius, sphereLen, Nx, Ntheta, varargin)

if (2*sphereLen > length)
    error('2*sphereLen must be less than the length');
end

quart = checkOptions({'Quarter'}, varargin);

if (mod(Nx, 2) ~= 0)
    error('The number of lengthwise panels (Nx) must be even');
end

if (quart)
    if (mod(Ntheta, 2) ~= 0)
        error('To generate a quarter geometry, the number of rotational panels (Ntheta) must be even');
    end
    dx = length/Nx;
    x = 0:dx:length/2;
    dtheta = pi/Ntheta;
    theta = -(pi/2):dtheta:0;
    
    Nx = Nx/2;
    Ntheta = Ntheta/2;
else
%     dx = length/Nx;
%     x = -length/2:dx:length/2;
    x = -length/2*cos(0:pi/Nx:pi);
    dtheta = pi/Ntheta;
    theta = -pi:dtheta:0;
end

ctheta = cos(theta);
stheta = sin(theta);


% make sure the corner of the cone lies on an x value.
[val ic] = min(abs(x - sphereLen));
sphereLen = x(ic);

slope = radius/sphereLen;


pans(Nx*Ntheta, 1) = Panel;

r = 0;


for n = 1:Nx
    if (x(n+1) >= (-length/2 + sphereLen) && (x(n+1) <= (length/2 - sphereLen)))
        rp1 = radius;
    else
        rp1 = radius*sqrt(1 - ((abs(x(n+1)) - (length/2 - sphereLen))/sphereLen)^2);
    end

    for m = 1:Ntheta
        verts = zeros(4,3);
        verts(1,:) = [x(n) r*ctheta(m) r*stheta(m)]; % top left
        verts(2,:) = [x(n) r*ctheta(m+1) r*stheta(m+1)]; % top right
        verts(3,:) = [x(n+1) rp1*ctheta(m+1) rp1*stheta(m+1)]; % bottom right
        verts(4,:) = [x(n+1) rp1*ctheta(m) rp1*stheta(m)]; % bottom left
            
        pans((n-1)*Ntheta + m) = Panel(verts);
    end
        r = rp1;
end

if (quart)
    geo = PanelGeo(pans, 'Xsym', 'Ysym');
else
    geo = PanelGeo(pans);
end
