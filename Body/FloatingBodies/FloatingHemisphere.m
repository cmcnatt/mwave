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
classdef FloatingHemisphere < FloatingBody
    
    properties (Access = protected)
        radius;
    end

    properties (Dependent)
        Radius;
    end
    
    methods
        % Constructor 
        function [fb] = FloatingHemisphere(rho, radius, varargin)
            % Floating half-submerged sphere
            %
            % rho - fluid density
            % radius - radius of the sphere
            % optional arguments - number or panels in theta (Ntheta), and
            % phi (Nphi)
            
            fb = fb@FloatingBody();
            fb.cg = [0 0 0];

            fb.m = FloatingHemisphere.MassMatrix(rho, radius);
            fb.position = [0 0 0];
            
            fb.radius = radius;
            
            if (~isempty(varargin))
                if (length(varargin) < 2)
                    error('There must be two optional inputs to FloatingHemisphere: Ntheta, Nphi');
                end
                Ntheta = varargin{1};
                Nphi = varargin{2};
                
                opts = checkOptions({{'UseSym'}, {'NoInt'}}, varargin);
                optsIn = {};
                n = 0;
                if (opts(1))
                    n = n + 1;
                    optsIn{n} = 'Quarter';
                end
                if (opts(2))
                    n = n + 1;
                    optsIn{n} = 'NoInt';
                    fb.iSurfPan = 0;
                else
                    fb.iSurfPan = 1;
                end
                
                panGeo = makePanel_hemisphere(fb.radius, Ntheta, Nphi, optsIn{:});
                fb.panelGeo = panGeo;
                fb.iLowHi = 0;
            end
            
            pts = CirContSurf.Compute(fb.cg, radius, 30);
            pts = [pts; radius 0 0];
            fb.wpSec = pts(:, 1:2);
        end
    end
    
    methods
        function [r] = get.Radius(fb)
            r = fb.radius;
        end
        
        function [] = MakePanelGeometry(fb, Ntheta, Nphi, varargin)
            opts = checkOptions({'UseSym'}, varargin);
            
            if (opts(1))
                panGeo = makePanel_hemisphere(fb.radius, Ntheta, Nphi, 'Quarter');
            else
                panGeo = makePanel_hemisphere(fb.radius, Ntheta, Nphi);
            end

            fb.panelGeo = panGeo;
            fb.iLowHi = 0;
        end
    end
    
    methods (Access = protected)
        function [] = onModifyCg(fb, cg)
            mass = fb.m(1);
            wetVol = 2/3*pi*fb.radius^3;
            rho = mass/wetVol;
            zcg = cg(3);

            fb.m = FloatingHemisphere.MassMatrix(rho, fb.radius, zcg);
        end
    end
    
    methods (Static)
        function [M] = MassMatrix(rho, radius, varargin)
            if (~isempty(varargin))
                % zcg is the vertical position in body coordinates, 
                % where z=0 is the calm water free surface
                zcg0 = 0;
                zcg = varargin{1};
                delzcg = zcg - zcg0;
            else
                delzcg = 0;
            end
            % rho is the fluid density
            % body mass equals buoyant force
            wetVol = 2/3*pi*radius^3;
            m = rho*wetVol;

            Izz = 2/5*m*radius^2;
            
            Ixx = Izz;
            Iyy = Izz;
            if (delzcg ~= 0)
                Ixx = Ixx + m*delzcg^2;
                Iyy = Iyy + m*delzcg^2;
            end

            M = zeros(6,6);

            M(1,1) = m;
            M(2,2) = m;
            M(3,3) = m;
            M(4,4) = Ixx;
            M(5,5) = Iyy;
            M(6,6) = Izz;
        end
    end
end