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
function [geo] = makePanel_spheroidEndHorCyl(len, radius, sphereLen, Nx, Ntheta, varargin)

if (2*sphereLen > len)
    error('2*sphereLen must be less than the length');
end

[opts, args] = checkOptions({{'Quarter'}, {'Hinge', 1}}, varargin);
quart = opts(1);

useHin = opts(2);
if (useHin)
    hingePos = args{2};
    
    Nhin = length(hingePos);
    dhin = 0.05*len;
    zhin = 0.3*radius;
    
    slophin = (radius - zhin)/dhin;
    hinxs = zeros(Nhin, 3);

    for m = 1:Nhin
        hinxs(m, 1) = hingePos(m) - dhin;
        hinxs(m, 2) = hingePos(m);
        hinxs(m, 3) = hingePos(m) + dhin;
    end
end

if (mod(Nx, 2) ~= 0)
    error('The number of lengthwise panels (Nx) must be even');
end

if (quart)
    % TODO: Quarter not implemented fully
    if (mod(Ntheta, 2) ~= 0)
        error('To generate a quarter geometry, the number of rotational panels (Ntheta) must be even');
    end
    dx = len/Nx;
    x = 0:dx:len/2;
    dtheta = pi/Ntheta;
    theta = -(pi/2):dtheta:0;
    
    Nx = Nx/2;
    Ntheta = Ntheta/2;
else
    %sphpct = sphereLen/len;
    %x1 = -len/2:len/Nx:len/2;
    x1 = -len/2*cos(0:pi/Nx:pi);
    
    if (~useHin)
        x = x1;
    else
        % make sure the corner of the cone lies on an x value.
        ihxs = zeros(size(hinxs));
        for m = 1:Nhin
            for n = 1:3
                [val, ihx] = min(abs(x1 - hinxs(m,n)));
                ihxs(m,n) = ihx;
            end
        end

        x = zeros(1,Nx+1);
        ihinx = zeros(1,Nx+1);

        for m = 1:3
            if (m == 1)
                starti = 1;
                startx = -len/2;
            else
                starti = ihxs(1,m-1)+1;
                startx = hinxs(1,m-1);
            end

            stopi = ihxs(1,m);
            nx = stopi - starti + 1;
            stopx = hinxs(1,m);    

            if (m == 1)
                del = pi/2/(nx-1);
                x(starti:stopi) = -(stopx - startx)*cos(0:del:pi/2) + stopx*ones(1,nx);
            else
                del = (stopx - startx)/nx;
                x(starti:stopi) = (startx + del):del:stopx;
                if (m == 2)
                    ihinx((starti-1):(stopi-1)) = ones(1,nx);
                else
                    ihinx((starti-1):(stopi-1)) = -ones(1,nx);
                end
            end  
        end

        for n = 2:Nhin    
            for m = 1:3
                if (m == 1)
                    starti = ihxs(n-1,3)+1;
                    startx = hinxs(n-1,3);
                else
                    starti = ihxs(n,m-1)+1;
                    startx = hinxs(n,m-1);
                end

                stopi = ihxs(n,m);
                nx = stopi - starti + 1;
                stopx = hinxs(n,m);      
                del = (stopx - startx)/nx;

                x(starti:stopi) = (startx + del):del:stopx;
                if (m == 2)
                    ihinx((starti-1):(stopi-1)) = ones(1,nx);
                elseif (m == 3)
                    ihinx((starti-1):(stopi-1)) = -ones(1,nx);
                end
            end
        end

        starti = ihxs(Nhin,3) + 1;
        stopi = Nx+1;
        nx = stopi - starti + 1;
        startx = hinxs(Nhin,3);
        stopx = len/2;    
        del = pi/2/nx;

        x(starti:stopi) = -(stopx - startx)*cos((pi/2+del):del:pi) + startx*ones(1,nx);
    end

    dtheta = pi/Ntheta;
    theta = -pi:dtheta:0;
end

ctheta = cos(theta);
stheta = sin(theta);


% make sure the corner of the cone lies on an x value.
[val, ic] = min(abs(x - sphereLen));
sphereLen = x(ic);

pans(Nx*Ntheta, 1) = Panel;

r = 0;

ihin = 0;
theta1 = theta;
hinz1 = -radius;

npan = 0;

for m = 1:Nx
    
    if (~useHin)
        ishin = false;
    else
        ishin = false;
        uphin = true;

        for n = 1:Nhin
            if ((((m + 1) > ihxs(n,1)) && ((m + 1) <= ihxs(n,3))) || (m >= ihxs(n,1)) && (m < ihxs(n,3)))
                ihin = n;
                ishin = true;
                if ((m + 1) > ihxs(n,2))
                    uphin = false;
                end
                break;
            end
        end
    end
        
    if (ishin)
        if (uphin)
            hinz2 = -radius + slophin*(x(m+1) - hinxs(ihin, 1));            
        else
            hinz2 = -zhin - slophin*(x(m+1) - hinxs(ihin, 2));
        end
        
        thetah22 = asin(hinz2/r);
        thetah21 = -pi - thetah22;
        
        % make sure the corner
        [~, iht21] = min(abs(theta - thetah21));
        [~, iht22] = min(abs(theta - thetah22));
                
        if (iht21 == iht22)
            theta2 = theta;
        else
            dthetah21 = -thetah22/(iht21-1);
            dthetah22 = (-thetah21 + thetah22)/(iht22 - iht21);
            dthetah23 = -thetah22/(Ntheta - iht22 + 1);

            theta2 = [-pi:dthetah21:thetah21, thetah21+dthetah22:dthetah22:thetah22, thetah22+dthetah23:dthetah23:0];
        end
        
        
        for n = 1:Ntheta
            verts = zeros(4,3);
            
            ztl = r*sin(theta1(n));
            ztr = r*sin(theta1(n+1));
            zbr = r*sin(theta2(n+1));
            zbl = r*sin(theta2(n));
            
            if (ztl < hinz1)
                ztl = hinz1;
            end
            
            if (ztr < hinz1)
                ztr = hinz1;
            end
            
            if (zbr < hinz2)
                zbr = hinz2;
            end
            
            if (zbl < hinz2)
                zbl = hinz2;
            end
            
            verts(1,:) = [x(m) r*cos(theta1(n)) ztl]; % top left
            verts(2,:) = [x(m) r*cos(theta1(n+1)) ztr]; % top right
            verts(3,:) = [x(m+1) r*cos(theta2(n+1)) zbr]; % bottom right
            verts(4,:) = [x(m+1) r*cos(theta2(n)) zbl]; % bottom left

            npan = npan + 1;
            pans(npan) = Panel(verts);
            if (pans(npan).Area == 0)
                
            end
            pans(npan).IsWet = true;
            pans(npan).IsBody = true;
            pans(npan).IsInterior = false;
        end
        
        hinz1 = hinz2;
        theta1 = theta2;
    else
        if (x(m+1) >= (-len/2 + sphereLen) && (x(m+1) <= (len/2 - sphereLen)))
            rp1 = radius;
        else
            rp1 = radius*sqrt(1 - ((abs(x(m+1)) - (len/2 - sphereLen))/sphereLen)^2);
        end

        for n = 1:Ntheta
            verts = zeros(4,3);
            verts(1,:) = [x(m) r*ctheta(n) r*stheta(n)]; % top left
            verts(2,:) = [x(m) r*ctheta(n+1) r*stheta(n+1)]; % top right
            verts(3,:) = [x(m+1) rp1*ctheta(n+1) rp1*stheta(n+1)]; % bottom right
            verts(4,:) = [x(m+1) rp1*ctheta(n) rp1*stheta(n)]; % bottom left

            npan = npan + 1;
            pans(npan) = Panel(verts);
            if (pans(npan).Area == 0)
                
            end
            pans(npan).IsWet = true;
            pans(npan).IsBody = true;
            pans(npan).IsInterior = false;
        end
        r = rp1;
        hinz1 = -radius;
        theta1 = theta;
    end
end

cnt = length(pans);
pansTop(cnt, 1) = Panel;
axis = [1 0 0];
angle = pi;

for n = 1:cnt
    pansTop(n) = Panel(pans(n).Vertices);
    pansTop(n).Rotate('AxisAngle', axis, angle);
    pansTop(n).IsBody = true;
    pansTop(n).IsWet = false;
    pansTop(n).IsInterior = false;
end

if (quart)
    geo = PanelGeo([pans; pansTop], 'Xsym', 'Ysym');
else
    geo = PanelGeo([pans; pansTop]);
end
