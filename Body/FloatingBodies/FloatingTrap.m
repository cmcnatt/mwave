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
        angF;
        angB;
        theta;
        trap;
        trapWet;
    end

    properties (Dependent)
        LengthTop;
        LengthBottom;
        Beam;
        Height;
        AngFront;
        AngBack;
        Draft;
        Theta;
    end
    
    methods
        % Constructor 
        function [fb] = FloatingTrap(rho, lenTop, lenBot, dLenF, beam, height, draft, varargin)
            fb = fb@FloatingBody();
            
            fb.position = [0 0 0];

            [m_, cg_, theta_, cb_] = FloatingTrap.MassMatrix(rho, lenTop, lenBot, dLenF, beam, height, draft);
            fb.m = m_;

            fb.cg = cg_;
            fb.cb = cb_;
            fb.theta = theta_;
            
            fb.lenTop = lenTop;
            fb.lenBot = lenBot;
            fb.beam = beam;
            fb.height = height;
            fb.draft = draft;
            
            if (~isempty(varargin))

                if (size(varargin,2) < 3)
                    error('There must be at least three optional inputs to FloatingFlap: Nx, Ny, Nz');
                end
                Nx = varargin{1};
                Ny = varargin{2};
                Nz = varargin{3};
                [panGeo, prof] =  makePanel_trap2(lenTop, lenBot, dLenF, beam, height, draft, Nx, Ny, Nz, varargin{:});
                
                opts = checkOptions({{'NoInt'}}, varargin);
                if (opts(1))
                    fb.iSurfPan = 0;
                else
                    fb.iSurfPan = 1;
                end
                fb.panelGeo = panGeo;
                fb.iLowHi = 0;
            end
            
            fb.angF = atan2(prof(4,1) - prof(1,1), (prof(4,2) - prof(1,2)));
            fb.angB = atan2(prof(2,1) - prof(3,1), (prof(2,2) - prof(3,2)));
            
            fb.trap = Trapazoid(prof(1:4,:));
            xZ = fb.trap.Xz(0);
            
            vwet = prof(1:4,:);
            vwet(1,:) = [xZ(1), 0];
            vwet(2,:) = [xZ(2), 0];
            
            fb.trapWet = Trapazoid(vwet);

            fb.wpSec = [xZ(2) -beam/2; xZ(2) beam/2; xZ(1) beam/2; xZ(1) -beam/2; xZ(2) -beam/2];
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
        
        function [a] = get.AngFront(fb)
            a = fb.angF;
        end
        
        function [a] = get.AngBack(fb)
            a = fb.angB;
        end
        
        function [thet] = get.Theta(fb)
            thet = fb.theta;
        end
                
        function [] = MakePanelGeometry(fb, Nx, Ny, Nz)
            warning('MakePanelGeometry not implemented');
        end
    end
    
    methods (Access = protected)
    end
        
    methods (Static)
        function [xZ] = Xz(lenTop, lenBot, dLenF, height, draft, z)
            
            mF = -dLenF/height;
            
            mB = (lenTop - lenBot - dLenF)/height;
            xzF = -lenTop/2-(height-draft-z)*mF;
            xzB = lenTop/2-(height-draft-z)*mB;
            
            xZ = [xzF xzB];
        end
        
        function [M, cg, theta, cb] = MassMatrix(rho, lenTop, lenBot, dLenF, beam, height, draft)
            % Area, cg, and moment of inertia of trapazoid assuming
            % constant density
            trp = Trapazoid.MakeTrap(lenTop, lenBot, dLenF, height, draft);
            Ag = trp.Area;
            Cg = trp.Centroid;
            Ig = trp.Interia;
            
            xZ = trp.Xz(0);
            lenZero = xZ(2) - xZ(1);
            verts = trp.Vertices;
            
            vwet = verts(1:4,:);
            vwet(1,:) = [xZ(1), 0];
            vwet(2,:) = [xZ(2), 0];
            
            trpWet = Trapazoid(vwet);
            
            % [Ag, Cg, Ig] = FloatingTrap.trapAreaCen(lenTop, lenBot, dLenF, height, draft);
           
            % Area, cb, and moment of inertia of trapzoid below water
            % assuming constant density
%             xZ = FloatingTrap.Xz(lenTop, lenBot, dLenF, height, draft, 0);
%             lenZero = xZ(2) - xZ(1);
%             dLenFZero = -lenTop/2+dLenF - xZ(1);
%            
%             [Ab, Cb, Ib] = FloatingTrap.trapAreaCen(lenZero, lenBot, dLenFZero, draft, draft);

            Ab = trpWet.Area;
            cb = trpWet.Centroid;
            
            % The cgx and cbx may not be the same, which will cause a
            % moment hydrostatic moment
            
            % Torque due to Cg and Cb offset
            Ts = (Cg(1)-cb(1))*Ab;
            
            % Moment arm of hydrostaic rotation
            Rarea = 1/3*((lenZero - Cg(1))^2 + (lenZero + Cg(1))^2);
            
            % Angle that the body would rotate to rest position.
            theta = Ts/Rarea;
            
            cg = [Cg(1) 0 Cg(2)];
            cb = [cb(1) 0 cb(2)];
            
            rhoBod = Ab/Ag*rho;
            m = rhoBod*Ag*beam;
                                    
            Iyy = rhoBod*beam*Ig;
            
            xCg = trp.Xz(cg(3));
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