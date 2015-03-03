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
function [geo, mass] = makePanelMass_doubleDuck(rho, r, lf, thetaf, lb, thetab, wid, Ntheta, Nr, Nwid, varargin)

opts = checkOptions({'Panel', 'Mass'}, varargin);

makePan = true;
makeMass = true;
if (opts)
    makePan = opts(1);
    makeMass = opts(2);
end

[xd, zd] = dduckProfile(r, lf, thetaf, lb, thetab, Ntheta);

% if the points are very close to the top, put them on the top to avoid
% having very small gaps/angles near the free surface
maxz = max([r, lf*cos(thetaf), lb*cos(thetab)]);

tol = 1e-4;
thetaTol = 10*pi/180; % 5 degrees

for n = 2:length(zd)
    if (abs(zd(n-1) - maxz) < tol)
        zd(n-1) = maxz;
        z1zero = true;
    else
        z1zero = false;
    end

    if (abs(zd(n) - maxz) < tol)
        zd(n) = maxz;
        z2zero = true;
    else
        z2zero = false;
    end
    
    if(xor(z1zero, z2zero))
        theta = atan2((zd(n) - zd(n-1)), (xd(n) - xd(n-1)));
        if (abs(theta) < thetaTol)
            if (z1zero)
                zd(n) = maxz;
            elseif (z2zero)
                zd(n-1) = maxz;
            end
        end
    end
end

for n = length(zd):-1:2
    if (abs(zd(n-1) - maxz) < tol)
        zd(n-1) = maxz;
        z1zero = true;
    else
        z1zero = false;
    end

    if (abs(zd(n) - maxz) < tol)
        zd(n) = maxz;
        z2zero = true;
    else
        z2zero = false;
    end
    
    if(xor(z1zero, z2zero))
        theta = atan2((zd(n) - zd(n-1)), (xd(n) - xd(n-1)));
        if (abs(theta) < thetaTol)
            if (z1zero)
                zd(n) = maxz;
            elseif (z2zero)
                zd(n-1) = maxz;
            end
        end
    end
end

r = 1:-1/Nr:0;

xds = r'*xd;
zds = r'*zd;

width = -wid/2:wid/Nwid:wid/2;


% left side (port side - front is facing the wave i.e. looking towards
% negative x)
pansL(Ntheta*Nr,1) = Panel;

% right side (starboard) 
pansR(Ntheta*Nr,1) = Panel;

% mass points
mps(Ntheta*Nr*Nwid,1) = MassPoint;

npan = 1;
nmas = 1;

for m = 1:Nr
    for n = 1:Ntheta
        verts = zeros(4,3);
        verts(1,:) = [xds(m,n) -wid/2 zds(m,n)];
        verts(2,:) = [xds(m+1,n) -wid/2 zds(m+1,n)];
        verts(3,:) = [xds(m+1,n+1) -wid/2 zds(m+1,n+1)];
        verts(4,:) = [xds(m,n+1) -wid/2 zds(m,n+1)];

        for o = 1:Nwid
            ypos1 = width(o)*ones(4,1);
            ypos2 = width(o+1)*ones(4,1);
            vertsM = [verts(:,1) ypos1 verts(:,3); verts(:,1) ypos2 verts(:,3)];

            mps(nmas) = QuadMassPoint(rho, vertsM);
            nmas = nmas + 1;
        end
        
        verts = flipud(verts);
        pansL(npan) = Panel(verts);
        pansL(npan).IsWet = true;
        pansL(npan).IsInterior = false;
        pansL(npan).IsBody = true;

        verts = zeros(4,3);
        verts(1,:) = [xds(m,n) wid/2 zds(m,n)];
        verts(2,:) = [xds(m+1,n) wid/2 zds(m+1,n)];
        verts(3,:) = [xds(m+1,n+1) wid/2 zds(m+1,n+1)];
        verts(4,:) = [xds(m,n+1) wid/2 zds(m,n+1)];


%        verts = flipud(verts);
        pansR(npan) = Panel(verts);
        pansR(npan).IsWet = true;
        pansR(npan).IsInterior = false;
        pansR(npan).IsBody = true;

        npan = npan + 1;
    end
end

pansS(Ntheta*Nwid, 1) = Panel;

npan = 1;
for m = 1:Nwid
    for n = 1:Ntheta
        verts = zeros(4,3);
        verts(1,:) = [xd(n) width(m) zd(n)];
        verts(2,:) = [xd(n) width(m+1) zd(n)];
        verts(3,:) = [xd(n+1) width(m+1) zd(n+1)];
        verts(4,:) = [xd(n+1) width(m) zd(n+1)];

        verts = flipud(verts);
        pansS(npan) = Panel(verts);
        pansS(npan).IsWet = true;
        pansS(npan).IsBody = true;
        if ((zd(n) == maxz) && (zd(n+1) == maxz))
            pansS(npan).IsInterior = true;
        else
            pansS(npan).IsInterior = false;
        end
        npan = npan + 1;
    end
end

geo = PanelGeo([pansL; pansS; pansR]);

mass = MassBody(mps);

end

function [x, z] = dduckProfile(r, lf, thetaf, lb, thetab, Ntheta)

nofront = false;
noback = false;

if (lf <= r)
    lf = r;
    thetaf = pi/2;
    nofront = true;
end

if (lb <= r)
    lb = r;
    thetab = pi/2;
    noback = true;
end

xc = 0;
zc = -r;

% the front point - based on th the front length and front angle.
xfp = -lf*sin(thetaf);
zfp = lf*cos(thetaf);

% the back point - based on the back length and the back angle
xbp = lb*sin(thetab);
zbp = lb*cos(thetab);

% tangent points on the circle from the front point on the top and the
% bottom
phif = acos(r/lf);
% xtfpt = -r*sin(thetaf - phif);
% ztfpt = r*cos(thetaf - phif);
% 
% xtfpb = -r*sin(thetaf + phif);
% ztfpb = r*cos(thetaf + phif);

% tangent points on the circle from the back point on the top and the
% bottom
phib = acos(r/lb);
% xtbpt = r*sin(thetab - phib);
% ztbpt = r*cos(thetab - phib);
% 
% xtbpb = r*sin(thetab + phib);
% ztbpb = r*cos(thetab + phib);

% create the top arc
dtheta = 2*pi/360;
theta1 = pi/2 + thetaf - phif;
theta2 = pi/2 - thetab + phib;

% no arc if they overlap
if (theta2 > theta1)
    xarct = [];
    zarct = [];
else
    arct = theta1:-dtheta:theta2;
    xarct = r*cos(arct);
    zarct = r*sin(arct);
end

% create the bottom arc
theta1 = pi/2 + thetaf + phif;
theta2 = pi/2 - thetab - phib + 2*pi;

% % no arc if they overlap
% if (theta2 < pi)
%     theta1 = theta1 - 2*pi;
% end

if (theta2 < theta1)
    xarcb = [];
    zarcb = [];
else
    arcb = theta2:-dtheta:theta1;
    xarcb = r*cos(arcb);
    zarcb = r*sin(arcb);
end


% assemble the points - from the front point clockwise
if (nofront && ~noback)
    xw = [xarct xbp xarcb];
    zw = [zarct zbp zarcb];
elseif (~nofront && noback)
    xw = [xfp xarct xarcb(2:end) xfp];
    zw = [zfp zarct zarcb(2:end) zfp];
elseif (nofront && noback)
    xw = [xarct xarcb(2:end)];
    zw = [zarct zarcb(2:end)];
else
    xw = [xfp xarct xbp xarcb xfp];
    zw = [zfp zarct zbp zarcb zfp];
end

[xw2, zw2] = resampleEven(xw, zw, Ntheta+1);

minDist = lb;
ibp = 1;
for n = 1:Ntheta
    dist = sqrt((xw2(n) - xbp)^2 + (zw2(n) - zbp)^2);
    if (dist < minDist)
        minDist = dist;
        ibp = n;
    end
end

xt = [xw2(1:ibp-1) xbp];
zt = [zw2(1:ibp-1) zbp];

xb = [xbp xw2(ibp+1:end)];
zb = [zbp zw2(ibp+1:end)];

[xt2, zt2] = resampleEven(xt, zt, length(xt));
[xb2, zb2] = resampleEven(xb, zb, length(xb));
   
x = [xt2 xb2(2:end)];
z = [zt2 zb2(2:end)];

end