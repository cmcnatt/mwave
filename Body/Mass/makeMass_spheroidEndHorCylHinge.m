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
function [cyl] = makeMass_spheroidEndHorCylHinge(rho, L, R, Ls, hingePos, Nx, Nr, Ntheta)

dx = L/Nx;
xl = -L/2+dx/2:dx:Nx*dx/2;

drc = R/Nr;
rc = drc/2:drc:(drc*Nr-drc/2);

dtheta = 2*pi/Ntheta;
sindth = sin(dtheta);
theta = 0:dtheta:(Ntheta-1)*dtheta;
sthet = sin(theta);
cthet = cos(theta);

Nhin = length(hingePos);
dhin = 0.05*L;
zhin = 0.3*R;

slophin = (R- zhin)/dhin;

hinxs = zeros(Nhin, 3);

for m = 1:Nhin
    hinxs(m, 1) = hingePos(m) - dhin;
    hinxs(m, 2) = hingePos(m);
    hinxs(m, 3) = hingePos(m) + dhin;
end

ihxs = zeros(size(hinxs));
for m = 1:Nhin
    for n = 1:3
        [~, ihx] = min(abs(xl - hinxs(m,n)));
        ihxs(m,n) = ihx;
    end
end

Ntot = Nx*Nr*Ntheta;

mps(Ntot,1) = MassPoint;

npoint = 1;

for l = 1:Nx
    x = xl(l);
    
    ishin = false;
    uphin = true;
    
    for n = 1:Nhin
        if (l > ihxs(n,1) && l < ihxs(n,3))
            ihin = n;
            ishin = true;
            if (l > ihxs(n,2))
                hinz = -zhin - slophin*(xl(l) - hinxs(ihin, 2));
            else
                hinz = -R + slophin*(xl(l) - hinxs(ihin, 1));
            end
            
            hR = -hinz./sthet;
            hr = zeros(Ntheta, Nr);
            hdr = zeros(Ntheta,1);
            for m = 1:Ntheta
                del = hR(m)/Ntheta;
                hdr(m) = del;
                hr(m,:) = del/2:del:(del*Nr - del/2);
            end
            
            thetah = zeros(1,4);
            thetah(1) = asin(-hinz/R);
            thetah(2) = pi - thetah(1);
            thetah(3) = pi + thetah(1);
            thetah(4) = 2*pi - thetah(1);
            break;
        end
    end
              
    if (abs(x) < (L/2 - Ls))
        r = rc;
        dr = drc;
    else
        xe = x - sign(x)*(L/2 - Ls);
        Re = R*sqrt(1 - (xe/Ls)^2);
        dre = Re/Nr;
        r = dre/2:dre:(dre*Nr-dre/2);
        dr = dre;
    end
    
    for m = 1:Nr
        for n = 1:Ntheta
            hinmass = false;
            if (ishin)
                if ((theta(n) > thetah(1) && (theta(n) < thetah(2))) || (theta(n) > thetah(3) && (theta(n) < thetah(4))))
                    hinmass = true;
                end
            end
            %hinmass = false;
            if (hinmass)
                % TODO: can't handle overlap of hinge with spheroidal end
                y = hr(n,m)*cthet(n);
                z = hr(n,m)*sthet(n);
                
                %vol = 2*(hr(n,m)*hdr(n) + slophin^2*dx^2)*sthet(n)^2*sindth*dx;
                vol = 2*hr(n,m)*hdr(n)*sthet(n)^2*sindth*dx;
                %volc = 1/2*((rc(m)+dr/2)^2 - (rc(m)-dr/2)^2)*dtheta*dx;
            else
                y = r(m)*cthet(n);
                z = r(m)*sthet(n);

                if (abs(x) < (L/2 - Ls))
                    vol = 1/2*((rc(m)+dr/2)^2 - (rc(m)-dr/2)^2)*dtheta*dx;
                else
                    x1 = xe + dx/2;
                    x2 = xe - dx/2;
                    rnd = 10^10;
                    val1 = sqrt(round(rnd*(1 - (x1/Ls)^2))/rnd);
                    val2 = sqrt(round(rnd*(1 - (x2/Ls)^2))/rnd);

                    vol = 1/2*Re*dtheta*dr*((x1/Ls*val1 - x2/Ls*val2) + real(asin(x1/Ls) - asin(x2/Ls)));

                    %vol = 1/2*dtheta*dr*dx*sqrt(Re^2 - (xe/Ls)^2);
                end
            end
            
            mps(npoint) = MassPoint(rho, vol, [x y z]);
            npoint = npoint + 1;
        end
    end
end

cyl = MassBody(mps);

end