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
function [box] = makeMass_box(rho, len, width, draft, freeboard, Nx, Ny, Nz)

dx = len/Nx;
xl = -len/2+dx/2:dx:Nx*dx/2;

dy = width/Ny;
yl = -width/2+dy/2:dy:Ny*dy/2;

hei = draft + freeboard;

dz = hei/Nz;
zl = -draft+dz/2:dz:Nz*dz/2;


vol = dx*dy*dz;
Ntot = Nx*Ny*Nz;

mps(Ntot,1) = MassPoint;

npoint = 1;
for l = 1:Nx
    for m = 1:Ny
        for n = 1:Nz
            x = xl(l);
            y = yl(m);
            z = zl(n);
            
            mps(npoint) = MassPoint(rho, vol, [x y z]);
            npoint = npoint + 1;
        end
    end
end

box = MassBody(mps);

end