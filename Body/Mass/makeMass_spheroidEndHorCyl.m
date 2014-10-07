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
function [cyl] = makeMass_spheroidEndHorCyl(rho, L, R, Ls, Nx, Nr, Ntheta)

dx = L/Nx;
xl = -L/2+dx/2:dx:Nx*dx/2;

drc = R/Nr;
rc = drc/2:drc:(drc*Nr-drc/2);

dtheta = 2*pi/Ntheta;
theta = 0:dtheta:(Ntheta-1)*dtheta;
sthet = sin(theta);
cthet = cos(theta);

Ntot = Nx*Nr*Ntheta;

mps = MassPoint.empty(Ntot, 0);

npoint = 1;

for l = 1:Nx
    x = xl(l);
    
    if (abs(x) < (L/2 - Ls))
        r = rc;
        dr = drc;
    else
        %xe = x - sign(x)*(L/2 - Ls);
        xe = (abs(x) - (L/2 - Ls))/Ls;
        dxe = dx/Ls;
        %Re = R*sqrt(1 - (xe/Ls)^2);
        Re = R*sqrt(1 - xe^2);
        dre = Re/Nr;
        r = dre/2:dre:(dre*Nr-dre/2);
        dr = dre;
        %dalpha2 = asin(mod(1/2*dx./r*R/Ls,1));
        alpha2 = asin((r+dr/2)./Re);
        alpha1 = asin((r-dr/2)./Re);
        dalpha = alpha2 - alpha1;
    end
    
    for m = 1:Nr
        for n = 1:Ntheta
            y = r(m)*cthet(n);
            z = r(m)*sthet(n);
            
            if (abs(x) < (L/2 - Ls))
                vol = 1/2*((rc(m)+dr/2)^2 - (rc(m)-dr/2)^2)*dtheta*dx;
            else
%                 x2 = xe + dxe/2;
%                 x1 = xe - dxe/2;
%                 rnd = 10^10;
%                 val2 = sqrt(round(rnd*(1 - (x2/Ls)^2))/rnd);
%                 val1 = sqrt(round(rnd*(1 - (x1/Ls)^2))/rnd);

%                 val2 = asin(round(rnd*x2)/rnd) + x2*sqrt(round(rnd*(1 - x2^2))/rnd);
%                 val1 = asin(round(rnd*x1)/rnd) + x1*sqrt(round(rnd*(1 - x1^2))/rnd);
                
                %vol = 1/2*R*dtheta*dr*(val2 - val1);
                
                %vol = R*Ls*dr*dtheta*dalpha(m) + r(m)*dr*dtheta*dx;
                vol = r(m)*dr*dtheta*dx; % this works ok - approximate as a cyliner
                %vol = (R^2/r(m) + r(m))*dr*dtheta*dx;
                
                
                %vol = 1/2*R*dtheta*dr*((x2/Ls*val1 - x1/Ls*val2) + real(asin(x2/Ls) - asin(x1/Ls)));
                %vol = 1/2*R*dtheta*dr*(abs(x2/Ls*val2 - x1/Ls*val1));
                
                %vol =2*r(m)*dr*dtheta*dx; % this one works
            end
            
            mps(npoint) = MassPoint(rho, vol, [x y z]);
            npoint = npoint + 1;
        end
    end
end

cyl = MassBody(mps);

end