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
function [geo] = makePanel_cylinder(radius, draft, height, Ntheta, Nr, Nz, varargin)

opts = checkOptions({{'Quarter'}, {'NoInt'}}, varargin);
quart = opts(1);
noInt = opts(2);

dr = radius/Nr;
r = radius:-dr:0;
dzn = draft/Nz;
Nzp = round((height - draft)/dzn);
if (height <= draft)
    z = 0:-dzn:-draft;
else
    dzp = (height - draft)/Nzp;
    z = [(height - draft):-dzp:0, -dzn:-dzn:-draft];
end
Nz = Nz + Nzp;
%z = 0:-dz:-draft;


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

if (z(1) == 0)    
    pans(Ntheta*(Nz + 2*Nr),1) = Panel;
    topZero = true;
else    
    if (noInt)
        pans(Ntheta*(Nz + 2*Nr),1) = Panel;
    else
        pans(Ntheta*(Nz + 3*Nr),1) = Panel;
    end
    topZero = false;
end

np = 0;

for n = 1:Ntheta
    
    % top
    for m = 1:Nr
        verts = zeros(4,3);
        verts(1,:) = [r(Nr-m+2)*ctheta(n) r(Nr-m+2)*stheta(n) z(1)];
        verts(2,:) = [r(Nr-m+1)*ctheta(n) r(Nr-m+1)*stheta(n) z(1)];
        verts(3,:) = [r(Nr-m+1)*ctheta(n+1) r(Nr-m+1)*stheta(n+1) z(1)];
        verts(4,:) = [r(Nr-m+2)*ctheta(n+1) r(Nr-m+2)*stheta(n+1) z(1)];

        np = np + 1;
        pans(np) = Panel(verts);
        if (topZero)
            pans(np).IsWet = true;
            if (noInt)
                pans(np).IsInterior = true;
            else
                pans(np).IsInterior = false;
            end
        else
            pans(np).IsWet = false;
        end
%         if (topZero)
%             pans(np).IsWet = true;
%         else
%             pans(np).IsWet = false;
%         end
%         if (topZero && ~noInt)
%             pans(np).IsInterior = true;
%         else
%             pans(np).IsInterior = false;
%         end
        pans(np).IsBody = true;
    end
    
    % interior for not top zero
    if (~noInt && ~topZero)
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
    
    % wall
    for m = 1:Nz
        verts = zeros(4,3);
        verts(1,:) = [cirx(n) ciry(n) z(m)];
        verts(2,:) = [cirx(n) ciry(n) z(m+1)];
        verts(3,:) = [cirx(n+1) ciry(n+1) z(m+1)];
        verts(4,:) = [cirx(n+1) ciry(n+1) z(m)];

        np = np + 1;
        pans(np) = Panel(verts);
        if (z(m) <= 0)
            pans(np).IsWet = true;
        else
            pans(np).IsWet = false;
        end
        pans(np).IsInterior = false;
        pans(np).IsBody = true;
    end

    % bottom
    for m = 1:Nr
        verts = zeros(4,3);
        verts(1,:) = [r(m)*ctheta(n) r(m)*stheta(n) -draft];
        verts(2,:) = [r(m+1)*ctheta(n) r(m+1)*stheta(n) -draft];
        verts(3,:) = [r(m+1)*ctheta(n+1) r(m+1)*stheta(n+1) -draft];
        verts(4,:) = [r(m)*ctheta(n+1) r(m)*stheta(n+1) -draft];

        np = np + 1;
        pans(np) = Panel(verts);
        pans(np).IsWet = true;
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