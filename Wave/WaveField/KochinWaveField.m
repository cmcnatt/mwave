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
classdef KochinWaveField < FuncWaveField
    % A wave field provides information (pressure, elevation, velocity,
    % significant wave height, spectrums) at discrete points in space
    % specified on either a grid (array) or at arbitrary (x,y,z) points.  
    %
    % Constructs a wave field of circular (cylindrical) waves radiating 
    % (or scattered) from a source based on the far-field Kochin function
    
    methods
        % Constructor
        function [wf] = KochinWaveField(rho, waves, isarray, varargin)
            if(checkOptions({'NoVel'}, varargin))
                compVel = false;
            else
                compVel = true;
            end
            if (isarray)
                args = cell(3,1);
                args(1:2) = varargin(1:2);
                args{3} = loc;
            else
                args = cell(2,1);
                args(1) = varargin(1);
                args{2} = loc;
            end
            wf = wf@FuncWaveField(rho, compVel, waves, @KochinWaveField.Compute, isarray, args{:});
        end
    end
    
    methods (Static)
       
        function [p, vel] = Compute(rho, compVel, waves, n, x, y, z, varargin)

            if (waves.IsPlane)
                loc = [0 0];
            else
                loc = waves.Origin;
            end
            
            kochinF = waves.KochinFunc(n);

            omega = waves.Omega(n);
            h = waves.H;
            k = waves.K(n);
            g = waves.G;

            x = x - loc(1);
            y = y - loc(2);

            r = sqrt(x.^2 + y.^2);
            theta = atan2(y,x);   
            
            [Nt, Mt] = size(theta);
            theta2 = theta;
            for n = 1:Nt
                for m = 1:Mt
                    if (theta2(n,m) < 0)
                        theta(n,m) = 2*pi+theta2(n,m);
                    end
                end
            end
            
            if (isinf(h))
                Kp = exp(k*z);

                Kpw = Kp;
            else
                Kp = cosh(k*(h + z))./cosh(k*h);

                Kpw = sinh(k*(h + z))./cosh(k*h);
            end

            kr = k*r;

            [f, df] = kochinF.Evaluate(theta);
            kreikr = (1./sqrt(kr)).*exp(-1i*kr);

            p = rho*g*Kp.*f.*kreikr;
            
            if (compVel)
                ur = g*k/omega*(1-1i./(2*kr)).*f.*kreikr;
                utheta = 1i*g/omega*Kp.*df.*kreikr;
                w = 1i*g*k/omega*Kpw.*p;

                cth = cos(theta);
                sth = sin(theta);
                u = ur.*cth - utheta.*sth;
                v = ur.*sth + utheta.*cth;

                vel = zeros([3, size(x)]);

                switch (ndims(vel))
                    case 2
                        vel(1,:) = u;
                        vel(2,:) = v;
                        vel(3,:) = w;
                    case 3
                        vel(1,:,:) = u;
                        vel(2,:,:) = v;
                        vel(3,:,:) = w;
                    case 4
                        vel(1,:,:,:) = u;
                        vel(2,:,:,:) = v;
                        vel(3,:,:,:) = w;
                    otherwise
                        error('velociy in plane wave has a unknown dimensions');
                end
            else
                vel = [];
            end
        end
    end
end