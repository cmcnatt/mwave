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
classdef FloatingCurvedFlap < FloatingBody
    
    properties (Access = protected)
        radTop;
        radBot;
        thick;
        draft;
    end

    properties (Dependent)
        RadiusTop;
        RadiusBot;
        Thickness;
        Draft;
    end
    
    methods
        % Constructor 
        function [fb] = FloatingCurvedFlap(rho, radTop, radBot, thick, draft, Nr, Ntheta, Nz, varargin)
            
            fb = fb@FloatingBody();

            fb.radTop = radTop;
            fb.radBot = radBot;
            fb.thick = thick;
            fb.draft = draft;
            
            fb.position = [0 0 0];
            fb.centRot = [0 0 -draft];

            fb.modes = ModesOfMotion;
            
            [fb.panelGeo, fb.wpSec] = makePanel_curvedFlap(radTop, radBot, thick, draft, Nr, Ntheta, Nz, varargin{:});
            fb.iLowHi = 0;
            
            [~, rhoBod, ~, ~, cg, cb, VolWet, VolBod] = computeHydroStatic(rho, fb.panelGeo, 0, fb.modes);
            fb.cg = cg;
            m = -rhoBod*VolBod;
            M = zeros(6,6);
            M(1,1) = m;
            M(2,2) = m;
            M(3,3) = m;
            fb.m = M;
            
            warning('Curved Flap mass computation not correct');
        end
        
        function [val] = get.RadiusTop(fb)
            val = fb.radTop;
        end
        
        function [val] = get.RadiusBot(fb)
            val = fb.radBot;
        end
        
        function [val] = get.Thickness(fb)
            val = fb.thick;
        end
        
        function [val] = get.Draft(fb)
            val = fb.draft;
        end
    end
end
