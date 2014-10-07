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
classdef CirContSurf < ControlSurface
    
    properties (Access = private)
        center;
        radius;
    end

    properties (Dependent)
        Center;
        Radius;
    end
    
    methods
        % Constructor
        function [surf] = CirContSurf(center, radius, npnts)
            cntr = zeros(1,3);
            cntr(1:2) = center(1:2);
            
            [pnts nrms ars] = CirContSurf.Compute(cntr, radius, npnts);
            
            surf = surf@ControlSurface(pnts, nrms, ars);            
            
            surf.center = cntr;
            surf.radius = radius;
        end
        
        % Center
        function [cntr] = get.Center(surf)
            cntr = surf.center;
        end
        
        % Radius
        function [rds] = get.Radius(surf)
            rds = surf.radius;
        end
    end
    
    methods (Static)
        function [pnts nrms ars] = Compute(center, radius, npnts)
            dTheta = 2*pi/npnts;
            ars = dTheta*radius*ones(npnts,1);

            theta = 0:dTheta:(2*pi-dTheta);

            pnts = zeros(npnts, 3);

            cTheta = cos(theta);
            sTheta = sin(theta);

            pnts(:,1) = center(1) + radius*cTheta;
            pnts(:,2) = center(2) + radius*sTheta;

            nrms = zeros(npnts, 3);

            nrms(:,1) = cTheta;
            nrms(:,2) = sTheta;
        end
    end
end