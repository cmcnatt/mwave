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
classdef IWaveField < handle
    % Wave field interface and abstract class.  
    
    properties (SetAccess = protected, GetAccess = protected)
        isarray;
        hasvel;
        rho;
        g;
        t;
        nT;
        h;
        points;
        x;
        y;
    end
    
    properties (SetAccess = private)
        IsArray;                    % Indicates whether the wave field is an array
        HasVelocity;                % Indicates whether the wave field containst velocity components
        Rho;                        % Water Density
        G;                          % Gravitational Constant
        T;                          % Periods
        H;                          % Depth
    end
    
    methods (Abstract)
        
        % Wave Field Values
        Pressure(wf, varargin)
        Elevation(wf, varargin)
        Velocity(wf, varargin)
        SigWaveHeight(wf, varargin)
        Spectra(wf, varargin)
        EnergyFlux(wf, surf, varargin)
        RemoveBodies(wf, bodies)
        
        % Overloaded Operators
        eq(wfa, wfb)
        ne(wfa, wfb)
        plus(wfa, wfb)
        uminus(wfin)
        minus(wfa, wfb)
        mtimes(a, b)
        times(a, b)
    end
    
    methods        
        % IsArray
        function [isAr] = get.IsArray(wf)
            isAr = wf.isarray;
        end
        
        % HasVelocity
        function [hv] = get.HasVelocity(wf)
            hv = wf.hasvel;
        end
        
        % Rho
        function [rh] = get.Rho(wf)
            rh = wf.rho;
        end
        
        % G
        function [g_] = get.G(wf)
            g_ = wf.g;
        end
        
        % T
        function [t_] = get.T(wf)
            t_ = wf.t;
        end
        
         % H
        function [h_] = get.H(wf)
            h_ = wf.h;
        end    
        
         % Locations where wave field evaluated.  Nx3 for field points.  [X Y Z] Meshgrid for field array.
        function [x_, y_] = FieldPoints(wf, varargin)
            if (wf.isarray)
                x_ = wf.x;
                y_ = wf.y;
            else
                surf = 0;
                for n = 1:length(varargin)
                    if (strcmp(varargin{n}, 'Surface'))
                        surf = 1;
                    end
                end
                
                if (surf)
                    cnt = 0;
                    if (wf.points(n,3) == 0)
                        cnt = cnt + 1;
                        indSurf(cnt) = n;
                    end
                    x_ = wf.points(indSurf);
                else
                    x_ = wf.points;
                end
            end
        end
        
        function [xo, ix, yo, iy] = FindClosestPoint(wave, xi, yi)
        % Finds the closest point and the indicies of that point in the
        % wave field
            xv = wave.x(1,:);
            yv = wave.y(:,1);

            if ((xi > max(xv)) || (yi > max(yv))  || (xi < min(xv)) || (yi < min(yv)))
                error('point is outside of the bounds of the wave field');
            end

            [buff, ix] = min(abs(xv - xi));
            [buff, iy] = min(abs(yv - yi));
            xo = wave.x(iy,ix);
            yo = wave.y(iy,ix);
        end
    end
end