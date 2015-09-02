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
classdef FloatingSphereEndCyl < FloatingBody
    
    properties (Access = protected)
        length;
        radius;
        sphereRad;
        hingePos;
        includeNotch;
    end

    properties (Dependent)
        Length;
        Radius;
        SphereRad;
        HingePos;
    end
    
    methods
        % Constructor 
        function [fb] = FloatingSphereEndCyl(rho, length, radius, sphereRad, hingePos, varargin)
            opts = checkOptions({{'Notch'}}, varargin);
            notch = opts(1);
            
            fb = fb@FloatingBody();
            fb.includeNotch = notch;
            
            if (2*sphereRad > length)
                error('2*sphereRad must be less than the length');
            end
            
            fb.length = length;
            fb.radius = radius;
            fb.sphereRad = sphereRad;
            
            fb.cg = [0 0 0];
            
            [modes, Nh] = FloatingSphereEndCyl.CreateHingeModes(length, hingePos);
            
            if (notch)
                fb.m = FloatingSphereEndCyl.MassMatrix(rho, length, radius, sphereRad, hingePos, 'Notch');
            else
                fb.m = FloatingSphereEndCyl.MassMatrix(rho, length, radius, sphereRad, hingePos);
            end
            fb.dpto = zeros(6 + Nh, 6 + Nh);
            fb.dpar = zeros(6 + Nh, 6 + Nh);
            fb.k = zeros(6 + Nh, 6 + Nh);

            fb.modes = modes;
            fb.nGen = Nh;
            fb.iGenMds = 28;
            fb.hingePos = hingePos;
            fb.writeFileMeth = @FloatingSphereEndCyl.writeHingePos;
            fb.writeParams = fb.hingePos;
            
            if (~isempty(varargin))
                if (size(varargin,2) < 2)
                    error('There must be two optional inputs to FloatingAttenuator: Nx, Ntheta');
                end
                Nx = varargin{1};
                Ntheta = varargin{2};
                if (notch)
                    fb.panelGeo = makePanel_spheroidEndHorCyl(fb.length, fb.radius, fb.sphereRad, Nx, Ntheta, 'Hinge', fb.hingePos);
                else
                    fb.panelGeo = makePanel_spheroidEndHorCyl(fb.length, fb.radius, fb.sphereRad, Nx, Ntheta);
                end
                fb.iLowHi = 0;
                
                if (Nh > 0)
                    fb.c = computeHydroStatic(rho, fb.panelGeo, 0, fb.modes);
                    fb.hasc = true;
                end
            end
                                    
            nx = 40;
            theta = 0:2*pi/nx:2*pi;
            x = -length/2*cos(theta);
            r = zeros(size(x));

            iends = (abs(x) > (length/2 - sphereRad));

            r(~iends) = radius*ones(1,sum(~iends));

            xend = abs(x) - (length/2 - sphereRad);
            rend = radius*sqrt(1 - xend(iends).^2./sphereRad^2);
            r(iends) = rend;

            posy = (theta <= pi);
            y = [r(posy), -r(~posy)];

            fb.wpSec = [x', y'];
        end
    end
    
    methods
        function [l] = get.Length(fb)
            l = fb.length;
        end
        
        function [r] = get.Radius(fb)
            r = fb.radius;
        end
        
        function [cpct] = get.SphereRad(fb)
            cpct = fb.conePct;
        end
        
        function [hp] = get.HingePos(fb)
            hp = fb.hingePos;
        end
        function [fb] = set.HingePos(fb, hp)
            [modes, Nh] = FloatingSphereEndCyl.CreateHingeModes(fb.length, hp);
            
            fb.m = FloatingSphereEndCyl.MassMatrix(fb.rho, fb.length, fb.radius, fb.sphereRad, hp);
            fb.d = zeros(6 + Nh, 6 + Nh);
            fb.k = zeros(6 + Nh, 6 + Nh);

            fb.modes = modes;
            fb.nGen = Nh;            
            
            fb.hingePos = hp;
            fb.writeParams = fb.hingePos;
        end
                
        function [] = MakePanelGeometry(fb, Nx, Ntheta)
            if (fb.includeNotch)
                fb.panelGeo = makePanel_spheroidEndHorCyl(fb.length, fb.radius, fb.sphereLen, fb.hingePos, Nx, Ntheta, 'Hinge', fb.hingePos);
            else
                fb.panelGeo = makePanel_spheroidEndHorCyl(fb.length, fb.radius, fb.sphereLen, fb.hingePos, Nx, Ntheta);
            end
            fb.iLowHi = 0;
        end
        
    end
    
    methods (Static)
        function [M] = MassMatrix(rho, length, radius, sphereRad, hingePos, varargin)
            notch = checkOptions({'Notch'}, varargin);
            Nx = 40;
            Nr = 20;
            Ntheta = 20;
            if (notch)
                cyl = makeMass_spheroidEndHorCylHinge(rho, length, radius, sphereRad, hingePos, Nx, Nr, Ntheta);
            else
                cyl = makeMass_spheroidEndHorCyl(rho, length, radius, sphereRad, Nx, Nr, Ntheta);
            end
                        
            modes = FloatingSphereEndCyl.CreateHingeModes(length, hingePos);
                        
            cyl.Modes = modes;
            
            M = cyl.MassMatrix;
        end
        
        function [modes, Nh] = CreateHingeModes(length, hingePos)
            if (~isempty(hingePos))
                [Nh, col] = size(hingePos);
                if (col > 1)
                    error('HingePos is an Nx1 vector of the x-position of the hinges');
                end
                for n = 1:Nh
                    if (abs(hingePos(n)) >= length/2)
                        error('The hinge position must be within the body length');
                    end
                end

                modes = ModesOfMotion();
                modes.Generalized = ones(1,Nh);

                for n = 1:Nh
                    hfunc = HingeYFunc();
                    hfunc.HingePos = [hingePos(n) 0];
                    hfuncs(n) = hfunc;
                end

                modes.GenMotFuncs = hfuncs;
            else
                modes = ModesOfMotion();
                Nh = 0;
            end
        end
    end
    
    methods (Static, Access = private)
        function [] = writeHingePos(path, fileName, params)
            hingePs = params{1};
            angle = params{3};
            
            filename = [path '\' fileName '_xhinge.dat'];
            fileID = fopen(filename, 'wt');

            Nhin = size(hingePs,1);
            fprintf(fileID, ['XHINGE.DAT file for hinge coordinates of hinged body used in ' fileName '\n']);
            fprintf(fileID, '0 %i \tISX,  NHIN\n', Nhin);
            for n = 1:Nhin
                fprintf(fileID, '%6.2f ', hingePs(n));
            end
            fprintf(fileID, '\tHINGPOS\n');
            fprintf(fileID, '%8.4f\tANGLE\n', pi/180*angle);

            fclose(fileID);
        end
    end
end