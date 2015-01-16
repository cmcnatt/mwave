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
classdef FloatingDuck < FloatingBody
    
    properties (Access = protected)
        rad;
        lenFront;
        angFront;
        lenBack;
        angBack;
        wid;
    end

    properties (Dependent)
        Radius;
    end
    
    methods
        % Constructor 
        function [fb] = FloatingDuck(rho, radius, lengthFront, angleFront, lengthBack, angleBack, width, varargin)
            % Salter type Duck, but can have 'wedges' on the front and/or back and at different angles
            %
            % rho - fluid density
            % radius - radius of the cylinder part
            % lengthFront - length from the center to tip of front 'wedge'
            % angleFront - angle from vertical toward negative x to 'spoke' out to front 'wedge'
            % lengthBack - length from the center to tip of back 'wedge'
            % angleBack - angle from vertical towards positive x to spoke out to back wedge
            % width - width of duck
            % optional arguments - number or panels in theta (Ntheta), r (Nr), and width (Nwid)
            
            fb = fb@FloatingBody();
            
            fb.cg = [0 0 0];
            fb.m = zeros(6,6);
            
            % extents
            maxz = max([radius, lengthFront*cos(angleFront), lengthBack*cos(angleBack)]);
            fb.position = [0 0 -maxz];
            
            fb.rad = radius;
            fb.lenFront = lengthFront;
            fb.angFront = angleFront;
            fb.lenBack = lengthBack;
            fb.angBack = angleBack;
            fb.wid = width;
            
            if (~isempty(varargin))
                if (length(varargin) ~= 3)
                    error('There must be three optional inputs to FloatingCylinder: Ntheta, Nr, Nz');
                end
                Ntheta = varargin{1};
                Nr = varargin{2};
                Nwid = varargin{3};
                [panGeo, mass] = makePanelMass_doubleDuck(rho, radius, lengthFront, angleFront, lengthBack, angleBack, width, Ntheta, Nr, Nwid);
                fb.panelGeo = panGeo;
                fb.iLowHi = 0;
                
                fb.m = mass.MassMatrix;
                fb.cg = mass.Cg;
                
                fb.iSurfPan = true;
            end
            
            maxxneg = -max([radius, lengthFront*sin(angleFront)]);
            maxxpos = max([radius, lengthFront*sin(angleFront)]);
            
            fb.wpSec = [maxxneg, -width/2; maxxneg, width/2; maxxpos width/2; maxxpos, -width/2; maxxneg, -width/2];
        end
    end
    
    methods
        function [r] = get.Radius(fb)
            r = fb.rad;
        end
                
        function [] = MakePanelGeometry(fb, Ntheta, Nr, Nwid)
            [panGeo, mass] = make_doubleDuck(fb.rad, fb.lenFront, fb.angFront, fb.lenBack, fb.angBack, fb.wid, Ntheta, Nr, Nwid);
            fb.panelGeo = panGeo;
            fb.iLowHi = 0;
            fb.m = mass.MassMatrix;
            fb.cg = mass.Cg;
            fb.iSurfPan = true;
        end
    end
end