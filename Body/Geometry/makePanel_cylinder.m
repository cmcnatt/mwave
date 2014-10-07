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
function [geo] = makePanel_cylinder(radius, draft, Ntheta, Nr, Nz, varargin)

quart = checkOptions({'Quarter'}, varargin);

dr = radius/Nr;
r = radius:-dr:0;
dz = draft/Nz;
z = 0:-dz:-draft;


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

% make outer surface
ctheta = cos(theta);
stheta = sin(theta);
cirx = radius*ctheta;
ciry = radius*stheta;

pans(Ntheta*(Nz + Nr),1) = Panel;

for n = 1:Ntheta
    for m = 1:Nz
        verts = zeros(4,3);
        verts(1,:) = [cirx(n) ciry(n) z(m)];
        verts(2,:) = [cirx(n) ciry(n) z(m+1)];
        verts(3,:) = [cirx(n+1) ciry(n+1) z(m+1)];
        verts(4,:) = [cirx(n+1) ciry(n+1) z(m)];

        pans((n-1)*(Nz+Nr) + m) = Panel(verts);
    end

    for m = 1:Nr
        verts = zeros(4,3);
        verts(1,:) = [r(m)*ctheta(n) r(m)*stheta(n) -draft];
        verts(2,:) = [r(m+1)*ctheta(n) r(m+1)*stheta(n) -draft];
        verts(3,:) = [r(m+1)*ctheta(n+1) r(m+1)*stheta(n+1) -draft];
        verts(4,:) = [r(m)*ctheta(n+1) r(m)*stheta(n+1) -draft];

        pans((n-1)*(Nz+Nr) + Nz + m) = Panel(verts);
    end
end

if (quart)
    geo = PanelGeo(pans, 'Xsym', 'Ysym');
else
    geo = PanelGeo(pans);
end


end