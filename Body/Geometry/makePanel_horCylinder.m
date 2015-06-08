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
function [geo] = makePanel_horCylinder(radius, beam, Nr, Ntheta, Ny, varargin)

[opts, args] = checkOptions({{'SphereEnd', 1}, {'FlareEnd', 1}}, varargin);

if (opts(1))
    sphLen = args{1};
else
    sphLen = 0;
end

if (opts(2))
    flareRad = args{2};
else
    flareRad = 0;
end

if (2*sphLen > beam)
    error('2*sphLen must be less than the beam');
end

if (flareRad > 0)
    if (flareRad < radius)
        error('Flare radius must be greater than radius');
    end
    
    ang = pi/4;
    
    delR = flareRad - radius;
    
    RR = delR/(1- cos(ang));
    
    startY = beam/2 - RR*sin(ang);
    if (startY < 0)
        error('Flare radius too big for beam');
    end
    
    dr = flareRad/Nr;
    r = flareRad:-dr:0;    
else
    dr = radius/Nr;
    r = radius:-dr:0;
    
end



dtheta = 2*pi/Ntheta;
theta = 0:dtheta:2*pi;

ctheta = cos(theta);
stheta = sin(theta);

dy = beam/Ny;
y = -beam/2:dy:beam/2;

rad = radius;
rp1 = radius;

if (sphLen > 0)
    Nr = 0;
    rad = 0;
    
    y = -beam/2*cos(0:pi/Ny:pi);
end

if (flareRad > 0)
    rad = flareRad;
end

pans(Ntheta*(2*Nr + Ny),1) = Panel;
np = 0;

for n = 1:Ntheta
    
    if (sphLen == 0)
        for m = 1:Nr
            verts = zeros(4,3);
            verts(1,:) = [r(Nr-m+2)*ctheta(n) y(1) r(Nr-m+2)*stheta(n)];
            verts(2,:) = [r(Nr-m+1)*ctheta(n) y(1) r(Nr-m+1)*stheta(n)];
            verts(3,:) = [r(Nr-m+1)*ctheta(n+1) y(1) r(Nr-m+1)*stheta(n+1)];
            verts(4,:) = [r(Nr-m+2)*ctheta(n+1) y(1) r(Nr-m+2)*stheta(n+1)];

            np = np + 1;
            pans(np) = Panel(verts);
            pans(np).IsWet = true;
            pans(np).IsInterior = false;
            pans(np).IsBody = true;
        end
    end
    
    for m = 1:Ny
        
        if (sphLen ~= 0)
            if (y(m+1) >= (-beam/2 + sphLen) && (y(m+1) <= (beam/2 - sphLen)))
                rp1 = radius;
            else
                rp1 = radius*sqrt(1 - ((abs(y(m+1)) - (beam/2 - sphLen))/sphLen)^2);
            end
        end
        
        if (flareRad ~= 0)
            if (y(m+1) >= -startY && (y(m+1) <= startY))
                rp1 = radius;
            else
                delY = abs(y(m+1)) - startY;
                thee = asin(delY/RR);
                rp1 = radius + RR*(1-cos(thee));
                %rp1 = radius+delY.^2;
            end            
        end
        
        verts = zeros(4,3);
        verts(1,:) = [rad*ctheta(n) y(m) rad*stheta(n)];
        verts(2,:) = [rp1*ctheta(n) y(m+1) rp1*stheta(n)];
        verts(3,:) = [rp1*ctheta(n+1) y(m+1) rp1*stheta(n+1)];
        verts(4,:) = [rad*ctheta(n+1) y(m) rad*stheta(n+1)];

        np = np + 1;
        pans(np) = Panel(verts);
        pans(np).IsWet = true;
        pans(np).IsInterior = false;
        pans(np).IsBody = true;
        
        rad = rp1;
    end

    if (sphLen == 0)
        for m = 1:Nr
            verts = zeros(4,3);
            verts(1,:) = [r(m)*ctheta(n) y(end) r(m)*stheta(n)];
            verts(2,:) = [r(m+1)*ctheta(n) y(end) r(m+1)*stheta(n)];
            verts(3,:) = [r(m+1)*ctheta(n+1) y(end) r(m+1)*stheta(n+1)];
            verts(4,:) = [r(m)*ctheta(n+1) y(end) r(m)*stheta(n+1)];

            np = np + 1;
            pans(np) = Panel(verts);
            pans(np).IsWet = true;
            pans(np).IsInterior = false;
            pans(np).IsBody = true;
        end
    end
end

geo = PanelGeo(pans);

end