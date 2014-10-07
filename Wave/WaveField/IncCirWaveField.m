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
classdef IncCirWaveField < FuncWaveField
    % A wave field provides information (pressure, elevation, velocity,
    % significant wave height, spectrums) at discrete points in space
    % specified on either a grid (array) or at arbitrary (x,y,z) points. 
    %
    % Constructs a wave field of incident waves (not plane waves) that are
    % specified by at bessel function summation
    %
    %   
    %   \-
    %   /_ f_m*J_m(k*r)e^(i*m*(theta - beta))
    %   
    % The summation is in general from negative infinity to infinity, but
    % is truncated to M (-M to M). The coefficients (f_m) represent the 
    % directional spread of the incident wave. They are either
    % entered directly or are computed from a function. 
    %
    % For plane incident waves, f_m = (-i)^m
    %
    % The wave component direction is the primary wave direction
    
    methods
        % Constructor
        function [wf] = IncCirWaveField(rho, waves, isarray, varargin)
            
            if(checkOptions({{'NoVel'}}, varargin))
                compVel = false;
            else
                compVel = true;
            end
            
            if (waves.IsPlane)
                error('Input waves must be CirWaves');
            end
            
            if (~waves.IsIncident)
                error('Input waves must be incident circular waves');
            end
            
            wf = wf@FuncWaveField(rho, compVel, waves, @IncCirWaveField.Compute, isarray, varargin{:});
        end
    end
    
    methods (Static)
       
            function [p, vel] = Compute(rho, compVel, waves, n, X, Y, Z, varargin)
                
                a = waves.A(n);
                h = waves.H(n);
                k = waves.K(n);
                beta = waves.Beta(n);
                epsilon = waves.Epsilon(n);
                g = waves.G;
                
                if (length(beta) < 1)
                    beta = 0;
                end
                
                f_m = waves.Coefs{n};
                orig = waves.Origin;
                
                x2 = X - orig(1);
                y2 = Y - orig(2);
                
                r = sqrt(x2.^2 + y2.^2);
                kr = k*r;
                theta = atan2(y2, x2);
                
                M = length(f_m);
                M = (M-1)/2;
                
                eta = f_m(M+1)*besselj(0, kr);
                for m = 1:M
                    Jm = besselj(m, kr);
                    cmthet = cos(m*(theta - beta));
                    smthet = sin(m*(theta - beta));
                    eta = eta + f_m(M+1+m)*Jm.*(cmthet + 1i*smthet) + f_m(M+1-m)*(-1)^m*Jm.*(cmthet - 1i*smthet);
                end
                
                eta = eta*a*exp(-1i*epsilon);
                
                if (isinf(h))
                    Kp = exp(k*Z);

                    Kpw = Kp;
                else
                    Kp = cosh(k*(h + Z))./cosh(k*h);

                    Kpw = sinh(k*(h + Z))./cosh(k*h);
                end

                if (n)
                    p = rho*g*Kp.*eta;
                else
                    p = -rho*g*Z + rho*g*Kp.*eta;
                end

                if (compVel)
                    % TODO: need to compute velocity
%                     u = g*kcbeta/omega*Kp.*eta;
%                     v = g*ksbeta/omega*Kp.*eta;
%                     w = g*k/omega*Kpw*1i.*eta;
                    u = 0;
                    v = 0;
                    w = 0;

                    vel = zeros([3, size(X)]);

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