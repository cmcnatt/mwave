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
function [geo, wpSec] = makePanel_curvedFlap(radTop, radBot, thick, draft, Nr, Ntheta, Nz, varargin)

[opts, args] = checkOptions({{'arc', 1}}, varargin);

theta0 = pi;
if opts(1)
    theta0 = pi/180*args{1};
end

dtheta = theta0/Ntheta;
theta = -theta0/2:dtheta:theta0/2;
dz = draft/Nz;
z = 0:-dz:-draft;
dr = thick/Nr;

slope = (radBot - radTop)/draft;

front = zeros(length(z), length(theta), 2);
back = zeros(length(z), length(theta), 2);
sideR = zeros(length(z), Nr+1, 2);
sideL = zeros(length(z), Nr+1, 2);
top = zeros(length(theta), Nr+1, 2);
bottom = zeros(length(theta), Nr+1, 2);

for m = 1:length(z)
    rad = radTop - slope*z(m);
    
    radF = rad - thick;
    radB = rad;
    
    if (radF < 0)
        error('CurvedFlap radius too small for thickness');
    end
    
    front(m,:,1) = radF*cos(theta);
    front(m,:,2) = radF*sin(theta);
    back(m,:,1) = radB*cos(theta);
    back(m,:,2) = radB*sin(theta);
    
    rads = radF:dr:radB;
    
    sideR(m, :, 1) = rads*cos(theta(1));
    sideR(m, :, 2) = rads*sin(theta(1));
    sideL(m, :, 1) = rads*cos(theta(end));
    sideL(m, :, 2) = rads*sin(theta(end));
    
    if (1 == m) || (length(z) == m)
        surf = zeros(length(theta), length(rads), 2);
        for n = 1:length(theta)
            surf(n, :, 1) = rads*cos(theta(n));
            surf(n, :, 2) = rads*sin(theta(n));
        end
        if 1 == m
            top = surf;
        else
            bottom = surf;
        end
    end     
end

wpSec = [squeeze(front(1,:,:)); flipud(squeeze(back(1,:,:))); squeeze(front(1,1,:))'];

pans(2*Ntheta*Nz + 2*Nr*Nz + Ntheta*Nr, 1) = Panel;

np = 0;

% front
for m = 1:Nz
    for n = 1:Ntheta
        verts = zeros(4,3);
        verts(1,:) = [front(m, n, 1) front(m, n, 2) z(m)]; % top right
        verts(2,:) = [front(m, n+1, 1) front(m, n+1, 2) z(m)]; % top left
        verts(3,:) = [front(m+1, n+1, 1) front(m+1, n+1, 2) z(m+1)]; % bottom left
        verts(4,:) = [front(m+1, n, 1) front(m+1, n, 2) z(m+1)]; % bottom right

        np = np + 1;
        pans(np) = Panel(verts);
        pans(np).IsWet = true;
        pans(np).IsBody = true;
        pans(np).IsInterior = false;
    end
end

%frontPanGeo = PanelGeo(pans);

% back
for m = 1:Nz
    for n = 1:Ntheta
        verts = zeros(4,3);
        verts(1,:) = [back(m, n+1, 1) back(m, n+1, 2) z(m)]; % top right
        verts(2,:) = [back(m, n, 1) back(m, n, 2) z(m)]; % top left
        verts(3,:) = [back(m+1, n, 1) back(m+1, n, 2) z(m+1)]; % bottom left
        verts(4,:) = [back(m+1, n+1, 1) back(m+1, n+1, 2) z(m+1)]; % bottom right

        np = np + 1;
        pans(np) = Panel(verts);
        pans(np).IsWet = true;
        pans(np).IsBody = true;
        pans(np).IsInterior = false;
    end
end

%fbPanGeo = PanelGeo(pans);

% sideR
for m = 1:Nz
    for n = 1:Nr
        verts = zeros(4,3);
        verts(1,:) = [sideR(m, n+1, 1) sideR(m, n+1, 2) z(m)]; % top right
        verts(2,:) = [sideR(m, n, 1) sideR(m, n, 2) z(m)]; % top left
        verts(3,:) = [sideR(m+1, n, 1) sideR(m+1, n, 2) z(m+1)]; % bottom left
        verts(4,:) = [sideR(m+1, n+1, 1) sideR(m+1, n+1, 2) z(m+1)]; % bottom right

        np = np + 1;
        pans(np) = Panel(verts);
        pans(np).IsWet = true;
        pans(np).IsBody = true;
        pans(np).IsInterior = false;
    end
end

% sideL
for m = 1:Nz
    for n = 1:Nr
        verts = zeros(4,3);
        verts(1,:) = [sideL(m, n, 1) sideL(m, n, 2) z(m)]; % top right
        verts(2,:) = [sideL(m, n+1, 1) sideL(m, n+1, 2) z(m)]; % top left
        verts(3,:) = [sideL(m+1, n+1, 1) sideL(m+1, n+1, 2) z(m+1)]; % bottom left
        verts(4,:) = [sideL(m+1, n, 1) sideL(m+1, n, 2) z(m+1)]; % bottom right

        np = np + 1;
        pans(np) = Panel(verts);
        pans(np).IsWet = true;
        pans(np).IsBody = true;
        pans(np).IsInterior = false;
    end
end

% top
for m = 1:Ntheta
    for n = 1:Nr
        verts = zeros(4,3);
        verts(1,:) = [top(m, n+1, 1) top(m, n+1, 2) z(1)]; % top right
        verts(2,:) = [top(m+1, n+1, 1) top(m+1, n+1, 2) z(1)]; % top left
        verts(3,:) = [top(m+1, n, 1) top(m+1, n, 2) z(1)]; % bottom left
        verts(4,:) = [top(m, n, 1) top(m, n, 2) z(1)]; % bottom right

        np = np + 1;
        pans(np) = Panel(verts);
        pans(np).IsWet = true;
        pans(np).IsBody = true;
        pans(np).IsInterior = true;
    end
end

% bottom
for m = 1:Ntheta
    for n = 1:Nr
        verts = zeros(4,3);
        verts(1,:) = [bottom(m, n, 1) bottom(m, n, 2) z(end)]; % top right
        verts(2,:) = [bottom(m+1, n, 1) bottom(m+1, n, 2) z(end)]; % top left
        verts(3,:) = [bottom(m+1, n+1, 1) bottom(m+1, n+1, 2) z(end)]; % bottom left
        verts(4,:) = [bottom(m, n+1, 1) bottom(m, n+1, 2) z(end)]; % bottom right

        np = np + 1;
        pans(np) = Panel(verts);
        pans(np).IsWet = true;
        pans(np).IsBody = true;
        pans(np).IsInterior = false;
    end
end

geo = PanelGeo(pans);
