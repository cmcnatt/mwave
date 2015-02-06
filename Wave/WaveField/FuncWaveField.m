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
classdef FuncWaveField < WaveField & handle
    % A wave field provides information (pressure, elevation, velocity,
    % significant wave height, spectrums) at discrete points in space
    % specified on either a grid (array) or at arbitrary (x,y,z) points.  
    %
    % A FuncWaveField creates a WaveField from a generic WaveField function
    

    methods
        
        % Constructor
        function [wf] = FuncWaveField(rho, compVel, waves, waveFunc, isarray, varargin)
                        
            if (length(waves) > 1)
                error('Only one IWaves object as an input to a FuncWaveField');
            end
            
            if (~isa(waves, 'IWaves'))
                error('waves must be of type IWaves');
            end
                        
            nT = waves.Count;
            T = zeros(nT, 1);

            optArgs = cell(1);
            
            if (isarray)
                X = varargin{1};
                Y = varargin{2};
                Z = zeros(size(X));
                [ny, nx] = size(X);
                p = zeros(nT, ny, nx);
                if (compVel)
                    vel = zeros(nT, 3, ny, nx);
                else
                    vel = [];
                end
                if (length(varargin) > 2)
                    optArgs = varargin(3:end);
                end
            else
                points = varargin{1};
                X = points(:,1);
                Y = points(:,2);
                Z = points(:,3);
                np = size(points, 1);
                p = zeros(nT, np);
                if (compVel)
                    vel = zeros(nT, 3, np);
                else
                    vel = [];
                end
                if (length(varargin) > 1)
                    optArgs = varargin(2:end);
                end
            end
            
            for n = 1:nT
                T(n) = waves.T(n);
                for m = 1:n-1
                    if (waves.T(n) == T(m))
                        error('Repeated frequencies (periods) no allowed in wave component');
                    end
                end
                
                if(compVel)
                    [pnm, velnm] = waveFunc(rho, true, waves, n, X, Y, Z, optArgs{:});
                    if(isarray)
                        p(n, :, :) = pnm;
                        vel(n, :, :, :) = velnm;
                    else
                        p(n, :) = pnm;
                        vel(n, :, :) = velnm;                        
                    end
                else
                    pnm = waveFunc(rho, false, waves, n, X, Y, Z, optArgs{:});
                    if(isarray)
                        p(n, :, :) = pnm;
                    else
                        p(n, :) = pnm;                  
                    end
                end
            end
                        
            g = waves.G;
            h = waves.H;
            
            wf = wf@WaveField(rho, g, h, T, p, vel, isarray, varargin{:});
        end
    end
end