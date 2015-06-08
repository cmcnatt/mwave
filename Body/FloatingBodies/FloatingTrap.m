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
classdef FloatingTrap < FloatingBody
    
    properties (Access = protected)
        lenTop;
        lenBot;
        beam;
        height;
        draft;
        theta;
        surfArea;
    end

    properties (Dependent)
        LengthTop;
        LengthBottom;
        Beam;
        Height;
        Draft;
        Theta;
        SurfaceArea;
    end
    
    methods
        % Constructor 
        function [fb] = FloatingTrap(rho, lenTop, lenBot, dLenF, beam, height, draft, varargin)
            fb = fb@FloatingBody();
            
            [m_, cg_, theta_] = FloatingTrap.MassMatrix(rho, lenTop, lenBot, dLenF, beam, height, draft);
            fb.m = m_;
            fb.position = [0 0 cg_(3)];
            fb.cg = [cg_(1) 0 0];
            fb.theta = theta_;
            
            fb.lenTop = lenTop;
            fb.lenBot = lenBot;
            fb.beam = beam;
            fb.height = height;
            fb.draft = draft;
            
            if (~isempty(varargin))
                if (size(varargin,2) ~= 3)
                    error('There must be three optional inputs to FloatingFlap: Nx, Ny, Nz');
                end
                Nx = varargin{1};
                Ny = varargin{2};
                Nz = varargin{3};
                %panGeo = makePanel_box(len, beam, draft, Nx, Ny, Nz, 'Quarter');
                panGeo =  makePanel_trap(lenTop, lenBot, dLenF, beam, height, draft, Nx, Ny, Nz);
                panGeo.Translate([0 0 -cg_(3)]);
                fb.panelGeo = panGeo;
                fb.iLowHi = 0;
            end
            
            xZ = FloatingTrap.Xz(lenTop, lenBot, dLenF, height, draft, 0);
            fb.wpSec = [xZ(2) -beam/2; xZ(2) beam/2; xZ(1) beam/2; xZ(1) -beam/2; xZ(2) -beam/2];
            
            ang1 = atan2(height, dLenF);
            lenSide1 = abs(height/sin(ang1));
            
            ang2 = atan2(height, lenBot + dLenF - lenTop);
            lenSide2 = abs(height/sin(ang1));
            
            sa = (lenTop + lenBot + lenSide1 + lenSide2)*beam + (lenTop + lenBot)*height;
            
            fb.surfArea = sa;
        end
    end
    
    methods
        function [l] = get.LengthTop(fb)
            l = fb.lenTop;
        end
        
        function [l] = get.LengthBottom(fb)
            l = fb.lenBot;
        end
        
        function [b] = get.Beam(fb)
            b = fb.beam;
        end
        
        function [h] = get.Height(fb)
            h = fb.height;
        end
        
        function [d] = get.Draft(fb)
            d = fb.draft;
        end
        
        function [thet] = get.Theta(fb)
            thet = fb.theta;
        end
        
        function [sa] = get.SurfaceArea(fb)
            sa = fb.surfArea;
        end
        
        function [] = MakePanelGeometry(fb, Nx, Ny, Nz)
            panGeo = makePanel_box(fb.length, fb.beam, fb.draft, Nx, Ny, Nz, 'Quarter');
            panGeo.Translate(-fb.position);
            fb.panelGeo = panGeo;
            fb.iLowHi = 0;
        end
    end
    
    methods (Static)
        function [xZ] = Xz(lenTop, lenBot, dLenF, height, draft, z)
            
            mF = -dLenF/height;
            
            mB = (lenTop - lenBot - dLenF)/height;
            xzF = -lenTop/2-(height-draft-z)*mF;
            xzB = lenTop/2-(height-draft-z)*mB;
            
            xZ = [xzF xzB];
        end
        
        function [M, cg, theta] = MassMatrix(rho, lenTop, lenBot, dLenF, beam, height, draft)
            
            % Area, cg, and moment of inertia of trapazoid assuming
            % constant density
            [Ag, Cg, Ig] = FloatingTrap.trapAreaCen(lenTop, lenBot, dLenF, height, draft);
           
            % Area, cb, and moment of inertia of trapzoid below water
            % assuming constant density
            xZ = FloatingTrap.Xz(lenTop, lenBot, dLenF, height, draft, 0);
            lenZero = xZ(2) - xZ(1);
            dLenFZero = -lenTop/2+dLenF - xZ(1);
           
            [Ab, Cb, Ib] = FloatingTrap.trapAreaCen(lenZero, lenBot, dLenFZero, draft, draft);
            
            % The cgx and cbx may not be the same, which will cause a
            % moment hydrostatic moment
            
            % Torque due to Cg and Cb offset
            Ts = (Cg(1)-Cb(1))*Ab;
            
            % Moment arm of hydrostaic rotation
            Rarea = 1/3*((lenZero - Cg(1))^2 + (lenZero + Cg(1))^2);
            
            % Angle that the body would rotate to rest position.
            theta = Ts/Rarea;
            
            cg = [Cg(1) 0 Cg(2)];
            
            rhoBod = Ab/Ag*rho;
            m = rhoBod*Ag*beam;
                                    
            Iyy = rhoBod*beam*Ig;
            
            xCg = FloatingTrap.Xz(lenTop, lenBot, dLenF, height, draft, cg(3));
            lenCg = xCg(2) - xCg(1);
            
            Ixx = 1/12*m*(height^2+beam^2);
            Izz = 1/12*m*(lenCg^2+beam^2);
            
            M = zeros(6,6);

            M(1,1) = m;
            M(2,2) = m;
            M(3,3) = m;
            M(4,4) = Ixx;
            M(5,5) = Iyy;
            M(6,6) = Izz;
        end
    end
    
    methods (Static, Access = private)
        function [A, C, I] = trapAreaCen(lenTop, lenBot, dLenF, height, draft)
            p(1,:) = [-lenTop/2, height-draft];
            p(2,:) = [lenTop/2, height-draft];
            p(3,:) = [-lenTop/2+dLenF+lenBot, -draft];
            p(4,:) = [-lenTop/2+dLenF, -draft];
           
            h = height; 
            % tri 1 topLen as base
            b = lenTop;
            A1 = 1/2*b*h;
            C1 = 1/3*(p(1,:) + p(2,:) + p(3,:));
            I1 = 1/36*(b^3*h-b^2*h*A1+b*h*A1^2+b*h^3);
            
            % tri 2 botLen as base
            b = lenBot;
            A2 = 1/2*b*h;
            C2 = 1/3*(p(3,:) + p(4,:) + p(1,:));
            I2 = 1/36*(b^3*h-b^2*h*A2+b*h*A2^2+b*h^3);
            
            A = A1 + A2;
            
            C = 1/(A1 + A2)*(C1*A1 + C2*A2);
                        
            d1 = C1 - C;
            d1 = sqrt(d1(1)^2 + d1(2)^2);
            
            d2 = C2 - C;
            d2 = sqrt(d2(1)^2 + d2(2)^2);
            
            I1 = I1 + A1*d1^2;
            I2 = I2 + A2*d2^2;
            
            I = I1 + I2;           
            
        end
    end
end