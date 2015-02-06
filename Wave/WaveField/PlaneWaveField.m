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
classdef PlaneWaveField < FuncWaveField
    % A wave field provides information (pressure, elevation, velocity,
    % significant wave height, spectrums) at discrete points in space
    % specified on either a grid (array) or at arbitrary (x,y,z) points. 
    %
    % Constructs a wave field of linear sinusoidal (plane) waves specified
    % by the wave components.  Plane waves can have multiple wave
    % frequencies but all wave components must have the same direction and 
    % there cannot have multiple wave components at the same frequency.
    
    properties (Access = private)
        beta;
    end
    
    methods
        % Constructor
        function [wf] = PlaneWaveField(rho, waves, isarray, varargin)  
            if(checkOptions({{'NoVel'}}, varargin))
                compVel = false;
            else
                compVel = true;
            end
            if (isarray)
                args = cell(2,1);
                args(1:2) = varargin(1:2);
            else
                args = cell(1,1);
                args(1) = varargin(1);
            end
            % compute the dynamic pressure intead of the total pressure
            idynPres = cell(length(waves),1);
            for n = 1:length(waves)
                idynPres{n} = 1;
            end
            wf = wf@FuncWaveField(rho, compVel, waves, @PlaneWaveField.Compute, isarray, args{:});
            if (~waves.UniformDir)
                error('All the waves must have the same direction');
            end
            wf.beta = waves.Beta(1);
        end
    end
    
    methods (Static)
       
            function [p, vel] = Compute(rho, compVel, waves, n, x, y, z, varargin)
                
                a = waves.A(n);
                omega = waves.Omega(n);
                h = waves.H;
                k = waves.K(n);
                beta = waves.Beta(n);
                epsilon = waves.Epsilon(n);
                g = waves.G;
                
                if length(beta) < 1
                    beta = 0;
                end
                
                kcbeta = k*cos(beta);
                ksbeta = k*sin(beta);

                eta = a*exp(-1i*(kcbeta*x + ksbeta*y - epsilon*ones(size(x))));

                if (isinf(h))
                    Kp = exp(k*z);

                    Kpw = Kp;
                else
                    Kp = cosh(k*(h + z))./cosh(k*h);

                    Kpw = sinh(k*(h + z))./cosh(k*h);
                end

                p = rho*g*Kp.*eta;
                
%                 if (idynp)
%                     p = rho*g*Kp.*eta;
%                 else
%                     p = -rho*g*z + rho*g*Kp.*eta;
%                 end

                if (compVel)
                    u = g*kcbeta/omega*Kp.*eta;
                    v = g*ksbeta/omega*Kp.*eta;
                    w = g*k/omega*Kpw*1i.*eta;

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

