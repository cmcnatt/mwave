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
classdef ControlSurface
    
    properties (Access = private)
        points;
        norms;
        areas;
        is2dxsec;
    end

    properties (Dependent)
        Points;
        Norms;
        Areas;
        Is2dXSection;
    end
    
    methods
        
        % Constructor
        function [surf] = ControlSurface(pnts, nrms, ars)
            if (nargin ~= 0)
                [np ncp] = size(pnts);
                [nn ncn] = size(nrms);
                [na nca] = size(ars);
                
                if ((np ~= nn) || (np ~= na))
                    error('The number of points must equal the number of normals and areas');
                end
                
                if ((ncp ~= 3) || (ncn ~= 3) || (nca ~= 1))
                    error('Each point and normal must have 3 components and the area is a single value');
                end
                
                if (all(pnts(:,3) == 0) && all(nrms(:,3) == 0))
                    surf.is2dxsec = true;
                end
                
                r = sqrt(nrms(:,1).^2 + nrms(:,2).^2 + nrms(:,3).^2);
                
                for n = 1:length(r)
                    if (r(n) ~= 0)
                        if (abs(r(n) - 1) > 0.1)
                            error('All normals do not have unit length');
                        end
                    end
                end
                
                if (any(pnts(:,3) > 0))
                    error('For a WetSurface, all points must have z values less than or equal to zero');
                end
                
                surf.points = pnts;
                surf.norms = nrms;
                surf.areas = ars;
            end                
        end
        
        % Points
        function [pnts] = get.Points(surf)
            pnts = surf.points;
        end
        
        % Norms
        function [nrms] = get.Norms(surf)
            nrms = surf.norms;
        end
        
        % Areas
        function [ars] = get.Areas(surf)
            ars = surf.areas;
        end
        
        % Is2dXSection
        function [is2d] = get.Is2dXSection(surf)
            is2d = surf.is2dxsec;
        end
    end
end