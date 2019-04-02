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
function [geo] = makePanel_circlePlate(radius, depth, thick, Ntheta, Nr)

dr = radius/Nr;
r = radius:-dr:0;

if ~isempty(thick)
    zt = -depth + thick/2;
    zb = -depth - thick/2;
    pans(Ntheta*2*Nr,1) = Panel;
else
    zt = -depth;
    zb = [];
    pans(Ntheta*Nr,1) = Panel;
end

dtheta = 2*pi/Ntheta;
theta = 0:dtheta:2*pi;
ctheta = cos(theta);
stheta = sin(theta);

np = 0;

for n = 1:Ntheta
    
    % top
    for m = 1:Nr
        verts = zeros(4,3);
        verts(1,:) = [r(Nr-m+2)*ctheta(n) r(Nr-m+2)*stheta(n) zt];
        verts(2,:) = [r(Nr-m+1)*ctheta(n) r(Nr-m+1)*stheta(n) zt];
        verts(3,:) = [r(Nr-m+1)*ctheta(n+1) r(Nr-m+1)*stheta(n+1) zt];
        verts(4,:) = [r(Nr-m+2)*ctheta(n+1) r(Nr-m+2)*stheta(n+1) zt];

        np = np + 1;
        pans(np) = Panel(verts);
        pans(np).IsWet = true;
        pans(np).IsInterior = false;
        pans(np).IsBody = true;
    end

    % bottom
    if ~isempty(thick)
        for m = 1:Nr
            verts = zeros(4,3);
            verts(1,:) = [r(m)*ctheta(n) r(m)*stheta(n) zb];
            verts(2,:) = [r(m+1)*ctheta(n) r(m+1)*stheta(n) zb];
            verts(3,:) = [r(m+1)*ctheta(n+1) r(m+1)*stheta(n+1) zb];
            verts(4,:) = [r(m)*ctheta(n+1) r(m)*stheta(n+1) zb];

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