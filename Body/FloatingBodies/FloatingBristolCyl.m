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
classdef FloatingBristolCyl < FloatingBody
    
    properties (Access = protected)
        beam;
        radius;
        depth;
    end

    properties (Dependent)
        Beam;
        Radius;
        Depth;
    end
    
    methods
        % Constructor 
        function [fb] = FloatingBristolCyl(rho, radius, beam, depth, varargin)

            [opts, args] = checkOptions({{'SphereEnd', 1}, {'FlareEnd', 1}}, varargin);
            if (opts(1))
                sphLen = args{1};
            else
                sphLen = 0;
            end
            
            if (opts(2))
                flrRad = args{2};
            else
                flrRad = 0;
            end
            
            fb = fb@FloatingBody();

            fb.beam = beam;
            fb.radius = radius;
            if (depth < radius)
                error('For now the Bristol cylinder must be completely submerged - depth >= radius');
            end
            fb.depth = depth;
            
            fb.position = [0 0 -depth];
            
            fb.cg = [0 0 0];
            M = FloatingCylinder.MassMatrixCgCenter(rho, radius, beam, beam);
            I55 = M(6,6);
            I66 = M(5,5);
            M(6,6) = I66;
            M(5,5) = I55;
            fb.m = M;
            fb.dpto = zeros(6, 6);
            fb.dpar = zeros(6, 6);
            fb.k = zeros(6, 6);
            fb.kgm = zeros(6, 6);

            fb.modes = ModesOfMotion;

            
            if (~isempty(varargin))
                if (size(varargin,2) < 3)
                    error('There must be three optional inputs to FloatingBristolCyl: Nr, Ntheta, Ny');
                end
                Nr = varargin{1};
                Ntheta = varargin{2};
                Ny = varargin{3};
                
                if (sphLen > 0)
                    warning('Bristol Cylinder mass computation not correct');
                    fb.panelGeo = makePanel_horCylinder(radius, beam, Nr, Ntheta, Ny, 'SphereEnd', sphLen);
                elseif (flrRad > 0)
                    warning('Bristol Cylinder mass computation not correct');
                    fb.panelGeo = makePanel_horCylinder(radius, beam, Nr, Ntheta, Ny, 'FlareEnd', flrRad);
                else
                    fb.panelGeo = makePanel_horCylinder(radius, beam, Nr, Ntheta, Ny);
                end
                fb.iLowHi = 0;
            end
                                   

            fb.wpSec = [-radius -beam/2; -radius beam/2; radius beam/2; radius -beam/2; -radius -beam/2];
        end
        
        function [r] = get.Radius(fb)
            r = fb.radius;
        end
        
        function [h] = get.Beam(fb)
            h = fb.height;
        end
        
        function [d] = get.Depth(fb)
            d = fb.draft;
        end
    end
end
