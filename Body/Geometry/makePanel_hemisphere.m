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
function [geo] = makePanel_hemisphere(radius, Ntheta, Nphi, varargin)

opts = checkOptions({{'Quarter'}, {'NoInt'}}, varargin);
quart = opts(1);
noInt = opts(2);

dphi = pi/Nphi;
phi = -pi/2:dphi:pi/2;

Nr = ceil(1/dphi);
dr = radius/Nr;
r = radius:-dr:0;

if (quart)
    if (mod(Ntheta, 4) ~= 0)
        error('The number of panels in the theta direction must be divisible by 4 in order to generate a quarter geometry');
    end
    
    Ntheta = Ntheta/4;
    
    dtheta = pi/2/Ntheta;
    theta = 0:dtheta:(pi/2);
else
    dtheta = 2*pi/Ntheta;
    theta = 0:dtheta:2*pi;
end

ctheta = cos(theta);
stheta = sin(theta);
rcphi = radius*cos(phi);
rsphi = radius*sin(phi);


if (noInt)
    pans(Ntheta*Nphi, 1) = Panel;
else
    pans(Ntheta*(Nphi + Nr), 1) = Panel;
end

np = 0;

for n = 1:Ntheta
    
    % Interior
    if ~noInt
        for m = 1:Nr
            verts = zeros(4,3);
            verts(1,:) = [r(Nr-m+2)*ctheta(n) r(Nr-m+2)*stheta(n) 0];
            verts(2,:) = [r(Nr-m+1)*ctheta(n) r(Nr-m+1)*stheta(n) 0];
            verts(3,:) = [r(Nr-m+1)*ctheta(n+1) r(Nr-m+1)*stheta(n+1) 0];
            verts(4,:) = [r(Nr-m+2)*ctheta(n+1) r(Nr-m+2)*stheta(n+1) 0];

            np = np + 1;
            pans(np) = Panel(verts);
            pans(np).IsWet = true;
            pans(np).IsInterior = true;
            pans(np).IsBody = false;
        end
    end
    
    % outer surface
    for m = 1:Nphi
        verts = zeros(4,3);
        verts(1,:) = [rcphi(m)*ctheta(n) rcphi(m)*stheta(n) rsphi(m)];
        verts(2,:) = [rcphi(m)*ctheta(n+1) rcphi(m)*stheta(n+1) rsphi(m)];
        verts(3,:) = [rcphi(m+1)*ctheta(n+1) rcphi(m+1)*stheta(n+1) rsphi(m+1)];
        verts(4,:) = [rcphi(m+1)*ctheta(n) rcphi(m+1)*stheta(n) rsphi(m+1)];

        np = np + 1;
        pans(np) = Panel(verts);
        if rsphi(m+1) > 0
            pans(np).IsWet = false;
        else
            pans(np).IsWet = true;
        end
        pans(np).IsInterior = false;
        pans(np).IsBody = true;
    end
end

if (quart)
    geo = PanelGeo(pans, 'Xsym', 'Ysym');
else
    geo = PanelGeo(pans);
end


end