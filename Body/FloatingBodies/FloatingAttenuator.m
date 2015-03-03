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
classdef FloatingAttenuator < FloatingBody
    
    properties (Access = protected)
        length;
        radius;
        conePct;
    end

    properties (Dependent)
        Length;
        Radius;
        ConePercent;
    end
    
    methods
        % Constructor 
        function [fb] = FloatingAttenuator(rho, length, radius, conePct, varargin)
            fb = fb@FloatingBody();
            
            if (conePct > 0.25)
                error('Cone percentage cannot be greater that 0.25');
            end
            
            fb.length = length;
            fb.radius = radius;
            fb.conePct = conePct;
            
            fb.cg = [0 0 0];
            
            fb.m = FloatingAttenuator.MassMatrix(rho, length, radius, conePct);
            fb.d = zeros(7, 7);
            fb.dpto = zeros(7, 7);
            fb.dpar = zeros(7, 7);
            fb.k = zeros(7, 7);
            
            fb.modes = ModesOfMotion([1 1 1 1 1 1 1]);
            fb.nGen = 1;
            
            if (~isempty(varargin))
                if (size(varargin,2) ~= 2)
                    error('There must be two optional inputs to FloatingAttenuator: Nx, Ntheta');
                end
                Nx = varargin{1};
                Ntheta = varargin{2};
                fb.panelGeo = makePanel_attenuator(length, radius, conePct, Nx, Ntheta);
                fb.iLowHi = 0;
            end
            
            fb.iGenMds = 27;
            l = length/2;
            cl = conePct*length;
            r = radius;
            pts = [-l, 0; -l+cl, r; -cl, r; 0, 0; cl, r; l-cl, r; l, 0; l-cl, -r; cl, -r; 0, 0; -cl, -r; -l+cl, -r; -l, 0];
            fb.wpSec = pts;
        end
    end
    
    methods
        function [l] = get.Length(fb)
            l = fb.length;
        end
        
        function [r] = get.Radius(fb)
            r = fb.radius;
        end
        
        function [cpct] = get.ConePercent(fb)
            cpct = fb.conePct;
        end
                
        function [] = MakePanelGeometry(fb, Nx, Ntheta)
            fb.panelGeo = makePanel_attenuator(fb.length, fb.radius, fb.conePct, Nx, Ntheta, 'Quarter');
            fb.iLowHi = 0;
        end
    end
    
    methods (Static)
        function [M] = MassMatrix(rho, length, radius, conePct)
            % length of a single cylinder section
            lcyl = length/2 - 2*conePct*length;
            % 1/2 the cylinder is submerged
            Vcyl = 1/2*pi*radius^2*lcyl;
            % mass of a single cylinder
            mcyl = rho*Vcyl;
                        
            % moments for single cyl about its cg
            IxxCyl = mcyl*radius^2/2;
            IyyCyl = 1/12*mcyl*(3*radius^2 + lcyl^2);
            
            % moments for single cyl about origin
            d = length/4;
            IyyCyl = IyyCyl + mcyl*d^2;
            
            % length of a single cone
            lcon = conePct*length;
            % 1/2 of the cone is submerged
            Vcon = 1/2*1/3*pi*radius^2*lcon;
            mcon = rho*Vcon;
            
            % moments for single cone about its cg (cg is 1/4*height of
            % cone)
            IxxCon = 3/10*mcon*radius^2;
            IyyCon = 3/20*mcon*radius^2 + 3/80*mcon*lcon^2;
            
            % moment for single inner cone
            IyyConi = IyyCon + mcon*(3/4*lcon)^2;
            % moment for single outer cone
            IyyCono = IyyCon + mcon*(5/4*lcon + lcyl)^2; 
            
            % total mass
            m = 2*mcyl + 4*mcon;
            
            % total inertia
            Ixx = 2*IxxCyl + 4*IxxCon;
            Iyy = 2*IyyCyl + 2*IyyConi + 2*IyyCono;
            Izz = Iyy;
            % Generalized mode
            Igg = Iyy;

            M = zeros(7,7);

            % TODO: this is incorret - need couplying terms between heave
            % and flex
            M(1,1) = m;
            M(2,2) = m;
            M(3,3) = m;
            M(4,4) = Ixx;
            M(5,5) = Iyy;
            M(6,6) = Izz;
            M(7,7) = Igg;
        end
    end
end