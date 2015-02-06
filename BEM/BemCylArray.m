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
classdef BemCylArray
    % class to define a cylinder of field points for computation of the
    % diffraction transfer matrix
    
    properties (Access = private)
        rad;
        nTheta;
        nZ;
        cosSpacing;
    end

    properties (Dependent)
        Radius;         % The radius of the cylinder
        Ntheta;         % The number of points in the theta direction
        Nz;             % The number of points z direction
        CosSpacing;     % Indicates whether to use cosine spacing in z
    end
    
    methods
        % Constructor 
        function [array] = BemCylArray(radius, Ntheta, Nz)
            % Constructor
            if (nargin == 0)
                array.rad = 0;
                array.nTheta = 1;
                array.nZ = 1;
            else
                array.rad = radius;
                array.nTheta = Ntheta;
                array.nZ = Nz;
            end
            array.cosSpacing = true;
        end
        
        function [rad] = get.Radius(array)
            rad = array.rad;
        end
        function [array] = set.Radius(array, ra)
            if (ra >= 0)    
                array.rad = ra;
            else
                error('The radius must not be negative');
            end
        end
        
        function [nt] = get.Ntheta(array)
            nt = array.nTheta;
        end
        function [array] = set.Ntheta(array, nt)
            if (isInt(nt))
                array.nTheta = nt;
            else
                error('The number of theta points must be an integer');
            end
        end

        function [nz] = get.Nz(array)
            nz = array.nZ;
        end
        function [array] = set.Nz(array, nz)
            if (isInt(nz))
                array.nTheta = nz;
            else
                error('The number of z points must be an integer');
            end
        end
        
        function [cs] = get.CosSpacing(array)
            cs = array.cosSpacing;
        end
        function [array] = set.CosSpacing(array, cs)
            if (isBool(cs))
                array.cosSpacing = cs;
            else
                error('The CosSpacing must be a Bool');
            end
        end
        
        function [points] = GetPoints(array, h)
            nz = array.nZ;                       
            nth = array.nTheta;                   
            
            if (array.cosSpacing)
                z = -h*(1 - cos(0:pi/2/(nz-1):pi/2)); 
                z = round(z*10^6)/10^6;       
            else
                z = linspace(0, -h, nz);
            end
            
            points = makeCirWFPoints(array.rad, nth, z); 
        end
    end
end