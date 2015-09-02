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
function [geo, profile] = makePanel_trap2(lenTop, lenBot, dLenF, width, height, draft, Nx, Ny, Nz, varargin)


% points on panel must wrap counter-clockwise
% we'll start on the front face, top left corner - go across a row left to
% right, then down (just like reading)
% then the left right face, back face, then left face - circling around the
% box facing it - at each face,  go from left to right then down

[opts, args] = checkOptions({{'NzF', 1}, {'NzB', 1}, {'CosSpc'}, {'NoInt'}}, varargin);

if (opts(1))
    Nzf = args{1};
else
    Nzf = Nz;
end

if (opts(2))
    Nzb = args{2};
else
    Nzb = Nz;
end

if (opts(3))
    cosSpc = true;
else
    cosSpc = false;
end

noInt = opts(4);

if (cosSpc)
    xT = cosspace(-lenTop/2, lenTop/2, Nx+1);
    xB = cosspace(-lenTop/2+dLenF, -lenTop/2+dLenF+lenBot, Nx+1);
    
    y = cosspace(-width/2, width/2, Ny+1);
else
    dxT = lenTop/Nx;
    xT = -lenTop/2:dxT:lenTop/2;
    dxB = lenBot/Nx;
    xB = (-lenTop/2+dLenF):dxB:(-lenTop/2+dLenF+lenBot);

    dy = width/Ny;
    y = -width/2:dy:width/2;
end

% regular z
[z, xz, Nz2, izzero] = compZ(draft, height, xT, xB, Nz, cosSpc);

% z front
[zf, xzf, Nz2f] = compZ(draft, height, xT, xB, Nzf, cosSpc);

% z back
[zb, xzb, Nz2b] = compZ(draft, height, xT, xB, Nzb, cosSpc);

if (noInt)
    pans(Ny*(Nz2f + Nzb)  + 2*Nx*Nz2 + Ny*Nx,1) = Panel;
else
    pans(Ny*(Nz2f + Nzb)  + 2*Nx*Nz2 + Ny*Nx,1) = Panel;
end

np = 0;

profile = [xT(1) z(1); xT(end) z(1); xB(end) z(end); xB(1) z(end); xT(1) z(1)];

% back
for n = 1:Nz2b
    for m = 1:Ny
        verts = zeros(4,3);
 
        verts(1,:) = [xzb(end,n) y(m+1) zb(n)]; % top right
        verts(2,:) = [xzb(end,n) y(m) zb(n)]; % top left
        verts(3,:) = [xzb(end,n+1) y(m) zb(n+1)]; % bottom left
        verts(4,:) = [xzb(end,n+1) y(m+1) zb(n+1)]; % bottom right


        np = np + 1;
        pans(np) = Panel(verts);
        if (zb(n+1) < 0)
            pans(np).IsWet = true;
        else
            pans(np).IsWet = false;
        end
        pans(np).IsBody = true;
        pans(np).IsInterior = false;
    end
end

% right - going from +x to -x
for n = 1:Nz2
    for m = Nx:-1:1
        verts = zeros(4,3);

        verts(1,:) = [xz(m,n) y(end) z(n)]; % top right
        verts(2,:) = [xz(m+1,n) y(end) z(n)]; % top left
        verts(3,:) = [xz(m+1,n+1) y(end) z(n+1)]; % bottom left
        verts(4,:) = [xz(m,n+1) y(end) z(n+1)]; % bottom right


        np = np + 1;
        pans(np) = Panel(verts);
        if (z(n+1) < 0)
            pans(np).IsWet = true;
        else
            pans(np).IsWet = false;
        end
        pans(np).IsBody = true;
        pans(np).IsInterior = false;
    end
end

% front - going from +y to -y
for n = 1:Nz2f
    for m = Ny:-1:1
        verts = zeros(4,3);

        verts(1,:) = [xzf(1,n) y(m) zf(n)]; % top right
        verts(2,:) = [xzf(1,n) y(m+1) zf(n)]; % top left
        verts(3,:) = [xzf(1,n+1) y(m+1) zf(n+1)]; % bottom left
        verts(4,:) = [xzf(1,n+1) y(m) zf(n+1)]; % bottom right

        np = np + 1;
        pans(np) = Panel(verts);
        if (zf(n+1) < 0)
            pans(np).IsWet = true;
        else
            pans(np).IsWet = false;
        end
        pans(np).IsBody = true;
        pans(np).IsInterior = false;
    end
end

% left
for n = 1:Nz2
    for m = 1:Nx
        verts = zeros(4,3);

        verts(1,:) = [xz(m+1,n) y(1) z(n)]; % top right
        verts(2,:) = [xz(m,n) y(1) z(n)]; % top left
        verts(3,:) = [xz(m,n+1) y(1) z(n+1)]; % bottom left
        verts(4,:) = [xz(m+1,n+1) y(1) z(n+1)]; % bottom right
        
        np = np + 1;
        pans(np) = Panel(verts);
        if (z(n+1) < 0)
            pans(np).IsWet = true;
        else
            pans(np).IsWet = false;
        end
        pans(np).IsBody = true;
        pans(np).IsInterior = false;
    end
end

% bottom - looking up
for n = Nx:-1:1
    for m = 1:Ny
        verts = zeros(4,3);

        verts(1,:) = [xz(n+1,end) y(m+1) z(end)]; % top right
        verts(2,:) = [xz(n+1,end) y(m) z(end)]; % top left
        verts(3,:) = [xz(n,end) y(m) z(end)]; % bottom left
        verts(4,:) = [xz(n,end) y(m+1) z(end)]; % bottom right

        np = np + 1;
        pans(np) = Panel(verts);
        pans(np).IsWet = true;
        pans(np).IsBody = true;
        pans(np).IsInterior = false;
    end
end

% top - looking up
for n = Nx:-1:1
    for m = 1:Ny
        verts = zeros(4,3);

        verts(1,:) = [xz(n+1,1) y(m+1) z(1)]; % top right
        verts(2,:) = [xz(n,1) y(m+1) z(1)]; % bottom right
        verts(3,:) = [xz(n,1) y(m) z(1)]; % bottom left
        verts(4,:) = [xz(n+1,1) y(m) z(1)]; % top left

        np = np + 1;
        pans(np) = Panel(verts);
        pans(np).IsBody = true;
        if (z(1) == 0)
            pans(np).IsInterior = true;
            pans(np).IsWet = true;
        else
            pans(np).IsInterior = false;
            pans(np).IsWet = false;
        end
    end
end

% interior - looking up
if (~noInt)
    if (z(1) ~= 0)
        for n = Nx:-1:1
            for m = 1:Ny
                verts = zeros(4,3);

                verts(1,:) = [xz(n+1,izzero) y(m+1) 0]; % top right
                verts(2,:) = [xz(n,izzero) y(m+1) 0]; % bottom right
                verts(3,:) = [xz(n,izzero) y(m) 0]; % bottom left
                verts(4,:) = [xz(n+1,izzero) y(m) 0]; % top left

                np = np + 1;
                pans(np) = Panel(verts);
                pans(np).IsWet = true;
                pans(np).IsBody = false;
                pans(np).IsInterior = true;
            end
        end
    end
end

geo = PanelGeo(pans);

end

function [z, xz, Nz2, izzero] = compZ(draft, height, xT, xB, Nz, cosSpc)
    dz = draft/Nz;
    NzP = round((height-draft)/dz);
    
    if (cosSpc)
        zp = cosspace(height-draft,0,NzP+1);
        zn = cosspace(0,-draft,Nz+1);
    else
        zp = linspace(height-draft,0,NzP+1);
        zn = linspace(0,-draft,Nz+1);
    end
    z = [zp, zn(2:end)];
    Nz2 = Nz+NzP;
    izzero = NzP+1;

    ms = (xT-xB)./height;
    x0 = xT-(height-draft)*ms;

    for n = 1:length(z)
        xz(n,:) = x0+z(n)*ms;
    end
    xz = xz';
end