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
function [geo] = makePanel_box(length, width, draft, Nx, Ny, Nz, varargin)

quart = checkOptions({'Quarter'}, varargin);

if (quart)
    if ((mod(Nx, 2) ~= 0) || (mod(Ny, 2) ~= 0) )
        error('The number of panels in the x-direction and in the y-dirction must be divisible by 2 in order to generate a quarter geometry');
    end
    
    % In symmetry mode - all x and y values must be positive
    % start with the back main face - from y=0 to y=beam/2 at constant x
    % then do side face - from x=length/2 to x=0 at constant y
    % then do bottom face
    
    dx = length/Nx;
    x = length/2:-dx:0;
    dy = width/Ny;
    y = 0:dy:width/2;
    dz = draft/Nz;
    z = 0:-dz:-draft;
    
    Nx = Nx/2;
    Ny = Ny/2;
    
    pans(Ny*Nz + Nx*Nz + Ny*Nx,1) = Panel;
    
    % back
    for n = 1:Nz
        for m = 1:Ny
            verts = zeros(4,3);
            verts(1,:) = [x(1) y(m+1) z(n)]; % top right
            verts(2,:) = [x(1) y(m) z(n)]; % top left
            verts(3,:) = [x(1) y(m) z(n+1)]; % bottom left
            verts(4,:) = [x(1) y(m+1) z(n+1)]; % bottom right

            pans((n-1)*Ny + m) = Panel(verts);
        end
    end
    
    % left
    for n = 1:Nz
        for m = 1:Nx
            verts = zeros(4,3);
            verts(1,:) = [x(m+1) y(end) z(n)]; % top right
            verts(2,:) = [x(m) y(end) z(n)]; % top left
            verts(3,:) = [x(m) y(end) z(n+1)]; % bottom left
            verts(4,:) = [x(m+1) y(end) z(n+1)]; % bottom right

            pans(Ny*Nz + (n-1)*Nx + m) = Panel(verts);
        end
    end
    
    % bottom - looking up
    for n = 1:Nx
        for m = 1:Ny
            verts = zeros(4,3);
            verts(1,:) = [x(n) y(m+1) z(end)]; % top right
            verts(2,:) = [x(n) y(m) z(end)]; % top left
            verts(3,:) = [x(n+1) y(m) z(end)]; % bottom left
            verts(4,:) = [x(n+1) y(m+1) z(end)]; % bottom right

            pans(Ny*Nz + Nx*Nz + (n-1)*Ny + m) = Panel(verts);
        end
    end
    
    geo = PanelGeo(pans, 'Xsym', 'Ysym');
else
    % points on panel must wrap counter-clockwise
    % we'll start on the front face, top left corner - go across a row left to
    % right, then down (just like reading)
    % then the left right face, back face, then left face - circling around the
    % box facing it - at each face,  go from left to right then down

    dx = length/Nx;
    x = -length/2:dx:length/2;
    dy = width/Ny;
    y = -width/2:dy:width/2;
    dz = draft/Nz;
    z = 0:-dz:-draft;

    pans(2*Ny*Nz + 2*Nx*Nz + Ny*Nx,1) = Panel;

    % front
    for n = 1:Nz
        for m = 1:Ny
            verts = zeros(4,3);
            verts(1,:) = [x(end) y(m+1) z(n)]; % top right
            verts(2,:) = [x(end) y(m) z(n)]; % top left
            verts(3,:) = [x(end) y(m) z(n+1)]; % bottom left
            verts(4,:) = [x(end) y(m+1) z(n+1)]; % bottom right

            pans((n-1)*Ny + m) = Panel(verts);
        end
    end

    % right - going from +x to -x
    for n = 1:Nz
        for m = Nx:-1:1
            verts = zeros(4,3);
            verts(1,:) = [x(m) y(end) z(n)]; % top right
            verts(2,:) = [x(m+1) y(end) z(n)]; % top left
            verts(3,:) = [x(m+1) y(end) z(n+1)]; % bottom left
            verts(4,:) = [x(m) y(end) z(n+1)]; % bottom right

            pans(Ny*Nz + (n-1)*Nx + (Nx-m) + 1) = Panel(verts);
        end
    end

    % back - going from +y to -y
    for n = 1:Nz
        for m = Ny:-1:1
            verts = zeros(4,3);
            verts(1,:) = [x(1) y(m) z(n)]; % top right
            verts(2,:) = [x(1) y(m+1) z(n)]; % top left
            verts(3,:) = [x(1) y(m+1) z(n+1)]; % bottom left
            verts(4,:) = [x(1) y(m) z(n+1)]; % bottom right

            pans(Ny*Nz + Nx*Nz + (n-1)*Ny + (Ny-m) + 1) = Panel(verts);
        end
    end

    % left
    for n = 1:Nz
        for m = 1:Nx
            verts = zeros(4,3);
            verts(1,:) = [x(m+1) y(1) z(n)]; % top right
            verts(2,:) = [x(m) y(1) z(n)]; % top left
            verts(3,:) = [x(m) y(1) z(n+1)]; % bottom left
            verts(4,:) = [x(m+1) y(1) z(n+1)]; % bottom right

            pans(2*Ny*Nz + Nx*Nz + (n-1)*Nx + m) = Panel(verts);
        end
    end

    % bottom - looking up
    for n = Nx:-1:1
        for m = 1:Ny
            verts = zeros(4,3);
            verts(1,:) = [x(n+1) y(m+1) z(end)]; % top right
            verts(2,:) = [x(n+1) y(m) z(end)]; % top left
            verts(3,:) = [x(n) y(m) z(end)]; % bottom left
            verts(4,:) = [x(n) y(m+1) z(end)]; % bottom right

            pans(2*Ny*Nz + 2*Nx*Nz + (Nx-n)*Ny + m) = Panel(verts);
        end
    end

    geo = PanelGeo(pans);
end

end