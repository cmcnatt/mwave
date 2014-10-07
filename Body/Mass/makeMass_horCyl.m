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
function [cyl] = makeMass_horCyl(rho, L, R, Nx, Nr, Ntheta)

dx = L/Nx;
xl = -L/2+dx/2:dx:Nx*dx/2;

dr = R/Nr;
r = dr/2:dr:(dr*Nr-dr/2);

dtheta = 2*pi/Ntheta;
theta = 0:dtheta:(Ntheta-1)*dtheta;

Ntot = Nx*Nr*Ntheta;

mps(Ntot,1) = MassPoint;

npoint = 1;
for l = 1:Nx
    for m = 1:Nr
        for n = 1:Ntheta
            x = xl(l);
            y = r(m)*cos(theta(n));
            z = r(m)*sin(theta(n));
            
            vol = 1/2*((r(m)+dr/2)^2 - (r(m)-dr/2)^2)*dtheta*dx;
            
            mps(npoint) = MassPoint(rho, vol, [x y z]);
            npoint = npoint + 1;
        end
    end
end

cyl = MassBody(mps);

end