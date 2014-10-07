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
function [sphere] = makeMass_sphere(rho, R, Nr, Ntheta, Nphi)

dr = R/Nr;
r = dr/2:dr:(dr*Nr-dr/2);

dtheta = 2*pi/Ntheta;
theta = 0:dtheta:(Ntheta-1)*dtheta;

dphi = pi/Nphi;
phi = dphi/2:dphi:dphi*Nphi;

Ntot = Nr*Ntheta*Nphi;

mps(Ntot,1) = MassPoint;

npoint = 1;
for l = 1:Nr
    for m = 1:Ntheta
        for n = 1:Nphi
            x = r(l)*sin(phi(n))*cos(theta(m));
            y = r(l)*sin(phi(n))*sin(theta(m));
            z = r(l)*cos(phi(n));
            
            vol = 1/3*((r(l)+dr/2)^3 - (r(l)-dr/2)^3)*(2*sin(phi(n))*sin(dphi/2))*dtheta;
            
            mps(npoint) = MassPoint(rho, vol, [x y z]);
            npoint = npoint + 1;
        end
    end
end

sphere = MassBody(mps);

end