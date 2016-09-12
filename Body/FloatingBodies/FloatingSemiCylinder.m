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
classdef FloatingSemiCylinder < FloatingBody
    
    properties (Access = protected)
        rad;
        draft;
    end

    properties (Dependent)
        Radius;
        Draft;
    end
    
    methods
        % Constructor 
        function [fb] = FloatingSemiCylinder(rho, rad, draft, Nr, Ntheta, Nz, varargin)
            
            fb = fb@FloatingBody();

            fb.rad = rad;
            fb.draft = draft;
            
            fb.position = [0 0 0];
            fb.centRot = [0 0 -draft];

            fb.modes = ModesOfMotion;
            
            [fb.panelGeo] = makePanel_semiCylinder(rad, draft, Nr, Ntheta, Nz, varargin{:});
            fb.iLowHi = 0;
            
            [vol, cg] = computeVolume(fb.panelGeo);
            fb.cg = cg;
            m = -rho*vol;
            M = zeros(6,6);
            M(1,1) = m;
            M(2,2) = m;
            M(3,3) = m;
            fb.m = M;
            
            dTheta = pi/Ntheta;
            theta = [-pi/2:dTheta:pi/2, -pi/2];
           
            fb.wpSec = rad*[cos(theta).', sin(theta).'];
            
            warning('Semi Cylinder mass computation not correct');
        end
        
        function [val] = get.Radius(fb)
            val = fb.rad;
        end
        
        function [val] = get.Draft(fb)
            val = fb.draft;
        end
    end
end
