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
classdef FloatingFlap < FloatingBody
    
    properties (Access = protected)
        len;
        beam;
        draft;
    end

    properties (Dependent)
        Length;
        Beam;
        Draft;
    end
    
    methods
        % Constructor 
        function [fb] = FloatingFlap(rho, length, beam, draft, varargin)
            fb = fb@FloatingBody();
            % put the hinge (origin) at the bottom
            fb.position = [0 0 -draft];
            % put the cg in center of the panel.
            fb.cg = [0 0 draft/2];
            fb.m = FloatingFlap.MassMatrix(rho, length, beam, draft, fb.cg);
            
            len = length;
            fb.len = len;
            fb.beam = beam;
            fb.draft = draft;
            
            if (~isempty(varargin))
                if (size(varargin,2) ~= 3)
                    error('There must be three optional inputs to FloatingFlap: Nx, Ny, Nz');
                end
                Nx = varargin{1};
                Ny = varargin{2};
                Nz = varargin{3};
                %panGeo = makePanel_box(len, beam, draft, Nx, Ny, Nz, 'Quarter');
                panGeo = makePanel_box(len, beam, draft, Nx, Ny, Nz);
                panGeo.Translate(-fb.position);
                fb.panelGeo = panGeo;
                fb.iLowHi = 0;
            end
            
            fb.wpSec = [len/2 -beam/2; len/2 beam/2; -len/2 beam/2; -len/2 -beam/2; len/2 -beam/2];
        end
    end
    
    methods
        function [l] = get.Length(fb)
            l = fb.len;
        end
        
        function [b] = get.Beam(fb)
            b = fb.beam;
        end
        
        function [d] = get.Draft(fb)
            d = fb.draft;
        end
        
        function [] = MakePanelGeometry(fb, Nx, Ny, Nz)
            panGeo = makePanel_box(fb.length, fb.beam, fb.draft, Nx, Ny, Nz, 'Quarter');
            panGeo.Translate(-fb.position);
            fb.panelGeo = panGeo;
            fb.iLowHi = 0;
        end
    end
    
    methods (Static)
        function [M] = MassMatrix(rho, length, beam, draft, cg)
            V = length*beam*draft;

            m = rho*V;

            Ixx = 1/12*m*(beam^2 + draft^2);
            Iyy = 1/12*m*(length^2 + draft^2);
            Izz = 1/12*m*(length^2 + beam^2);
            
            if (cg(1) ~= 0)
                Iyy = Iyy + m*cg(1)^2;
                Izz = Izz + m*cg(1)^2;
            end
            
            if (cg(2) ~= 0)
                Ixx = Ixx + m*cg(2)^2;
                Izz = Izz + m*cg(2)^2;
            end
            
            if (cg(3) ~= 0)
                Ixx = Ixx + m*cg(3)^2;
                Iyy = Iyy + m*cg(3)^2;
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