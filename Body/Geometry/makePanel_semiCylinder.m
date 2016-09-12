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
function [geo] = makePanel_semiCylinder(rad, draft, Nr, Ntheta, Nz, varargin)

[opts, args] = checkOptions({{'flip'}}, varargin);

flip = opts(1);

dtheta = pi/Ntheta;
theta = -pi/2:dtheta:pi/2;

dz = draft/Nz;
z = 0:-dz:-draft;
dr = rad/Nr;
r = 0:dr:rad;

%pans(2*Nr*Nz + 2*Nr*Ntheta + Ntheta*Nr, 1) = Panel;

np = 0;

% front
for m = 1:Nz
    % 1
    for n = 1:Nr
        verts = zeros(4,3);
        verts(1,:) = [0, r(n), z(m)]; % top right
        verts(2,:) = [0, r(n+1), z(m)]; % top left
        verts(3,:) = [0, r(n+1), z(m+1)]; % bottom left
        verts(4,:) = [0, r(n), z(m+1)]; % bottom right

        np = np + 1;
        pans(np) = Panel(verts);
        pans(np).IsWet = true;
        pans(np).IsBody = true;
        pans(np).IsInterior = false;
    end
    
    % 2
    for n = 1:Nr
        verts = zeros(4,3);
        verts(1,:) = [0, -r(n+1), z(m)]; % top right
        verts(2,:) = [0, -r(n), z(m)]; % top left
        verts(3,:) = [0, -r(n), z(m+1)]; % bottom left
        verts(4,:) = [0, -r(n+1), z(m+1)]; % bottom right

        np = np + 1;
        pans(np) = Panel(verts);
        pans(np).IsWet = true;
        pans(np).IsBody = true;
        pans(np).IsInterior = false;
    end
end

% back
for m = 1:Nz
    % 1
    for n = 1:Ntheta
        verts = zeros(4,3);
        verts(1,:) = [rad*cos(theta(n+1)), rad*sin(theta(n+1)), z(m)]; % top right
        verts(2,:) = [rad*cos(theta(n)), rad*sin(theta(n)), z(m)]; % top left
        verts(3,:) = [rad*cos(theta(n)), rad*sin(theta(n)), z(m+1)]; % bottom left
        verts(4,:) = [rad*cos(theta(n+1)), rad*sin(theta(n+1)), z(m+1)]; % bottom right

        np = np + 1;
        pans(np) = Panel(verts);
        pans(np).IsWet = true;
        pans(np).IsBody = true;
        pans(np).IsInterior = false;
    end
end

% top
for m = 1:Nr
    for n = 1:Ntheta
        verts = zeros(4,3);
        verts(1,:) = [r(m+1)*cos(theta(n)), r(m+1)*sin(theta(n)), 0]; % top right
        verts(2,:) = [r(m+1)*cos(theta(n+1)), r(m+1)*sin(theta(n+1)), 0]; % top left
        verts(3,:) = [r(m)*cos(theta(n+1)), r(m)*sin(theta(n+1)), 0]; % bottom left
        verts(4,:) = [r(m)*cos(theta(n)), r(m)*sin(theta(n)), 0]; % bottom right

        np = np + 1;
        pans(np) = Panel(verts);
        pans(np).IsWet = true;
        pans(np).IsBody = true;
        pans(np).IsInterior = true;
    end
end

% bottom
for m = 1:Nr
    for n = 1:Ntheta
        verts = zeros(4,3);
        verts(1,:) = [r(m)*cos(theta(n)), r(m)*sin(theta(n)), -draft]; % top right
        verts(2,:) = [r(m)*cos(theta(n+1)), r(m)*sin(theta(n+1)), -draft]; % top left
        verts(3,:) = [r(m+1)*cos(theta(n+1)), r(m+1)*sin(theta(n+1)), -draft]; % bottom left
        verts(4,:) = [r(m+1)*cos(theta(n)), r(m+1)*sin(theta(n)), -draft]; % bottom right

        np = np + 1;
        pans(np) = Panel(verts);
        pans(np).IsWet = true;
        pans(np).IsBody = true;
        pans(np).IsInterior = false;
    end
end

geo = PanelGeo(pans);

if flip
    geo.Rotate([0 0 1], pi);
    geo.Translate([-rad, 0, 0]);
end

end