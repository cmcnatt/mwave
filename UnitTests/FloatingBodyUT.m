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
classdef FloatingBodyUT < matlab.unittest.TestCase

    methods (Test)
       
        function testCylinder(testCase)
            rho = 1000;
            r = 10;
            h = 20;
            d = 10;

            fb = FloatingCylinder(rho, r, h, d);
            
            testCase.verifyEqual(fb.Radius, r);
            testCase.verifyEqual(fb.Height, h);
            testCase.verifyEqual(fb.Draft, d);
            testCase.verifyEqual(fb.Cg, [0 0 0]);
            testCase.verifyEqual(fb.Zpos, -d+h/2);

            m = rho*pi*r^2*d;
            Izz = m*r^2/2;
            Ixx = 1/12*m*(3*r^2 + h^2);
            Iyy = Ixx;

            M = zeros(6,6);
            M(1,1) = m;
            M(2,2) = m;
            M(3,3) = m;
            M(4,4) = Ixx;
            M(5,5) = Iyy;
            M(6,6) = Izz;
            
            testCase.verifyEqual(fb.M, M, 'AbsTol', 1e-6);
            sect = fb.WaterPlaneSec;

            figure;
            subplot(2,1,1);
            plot(sect(:,1), sect(:,2));
            title({'FloatingBodyUT - testCylinder', ' ', 'Cylinder waterplane section'});
            axis equal;
            set(gca, 'xlim', 2*[-r r], 'ylim', 2*[-r r]);

            Ntheta = 12;
            Nr = 5;
            Nz = 8;

            fb.MakePanelGeometry(Ntheta, Nr, Nz);
            geo = fb.PanelGeo;

            subplot(2,1,2);
            plot(geo,'ShowSym');
            title('Cylinder panel geometry');
            axis equal;
        end
        
        function testFlap(testCase)
            rho = 1000;
            l = 2;
            b = 20;
            d = 10;

            Nx = 2;
            Ny = 20;
            Nz = 10;

            fb = FloatingFlap(rho, l, b, d, Nx, Ny, Nz);
            
            testCase.verifyEqual(fb.Length, l);
            testCase.verifyEqual(fb.Beam, b);
            testCase.verifyEqual(fb.Draft, d);
            testCase.verifyEqual(fb.Cg, [0 0 d/2]);
            testCase.verifyEqual(fb.Zpos, -d);

            m = rho*l*b*d;

            Ixx = 1/12*m*(b^2 + d^2) + m*(d/2)^2;
            Iyy = 1/12*m*(l^2 + d^2) + m*(d/2)^2;
            Izz = 1/12*m*(l^2 + b^2);

            M = zeros(6,6);

            M(1,1) = m;
            M(2,2) = m;
            M(3,3) = m;
            M(4,4) = Ixx;
            M(5,5) = Iyy;
            M(6,6) = Izz;
            
            testCase.verifyEqual(fb.M, M, 'AbsTol', 1e-6);

            sect = fb.WaterPlaneSec;

            figure;
            subplot(2,1,1);
            plot(sect(:,1), sect(:,2));
            title({'FloatingBodyUT - testFlap', ' ', 'Flap waterplane section'});
            axis equal;
            set(gca, 'xlim', [-b b], 'ylim', [-b b]);

            geo = fb.PanelGeo;

            subplot(2,1,2);
            plot(geo, 'ShowSym');
            title('Flap panel geometry');
            axis equal;
        end
        
        function testAttenuator(testCase)
            rho = 1000;
            l = 40;
            r = 2;
            conePct = 0.05;

            Nx = 40;
            Ntheta = 8;

            fb = FloatingAttenuator(rho, l, r, conePct, Nx, Ntheta);
            
            testCase.verifyEqual(fb.Length, l);
            testCase.verifyEqual(fb.Radius, r);
            testCase.verifyEqual(fb.ConePercent, conePct);
            testCase.verifyEqual(fb.Cg, [0 0 0]);
            testCase.verifyEqual(fb.Zpos, 0);

            testCase.verifyEqual(size(fb.M), [7 7]);
            testCase.verifyEqual(size(fb.Dpto), [7 7]);
            testCase.verifyEqual(size(fb.K), [7 7]);

            sect = fb.WaterPlaneSec;

            figure;
            subplot(2,1,1);
            plot(sect(:,1), sect(:,2));
            title({'FloatingBodyUT - testAttenuator', ' ', 'Attenuator waterplane section'});
            axis equal;
            set(gca, 'xlim', [-l l], 'ylim', [-l l]);

            geo = fb.PanelGeo;

            subplot(2,1,2)
            plot(geo,'ShowSym');
            title('Attenuator panel geometry');
            axis equal;
        end
        
        function testFloatSphereEndCylHinge(testCase)
            rho = 1000;
            l = 10;
            r = 0.5;

            sphereRad = 0.1*l;
            hingePos = [-l/6; l/6];

            Nx = 120;
            Ntheta = 16;

            fb = FloatingSphereEndCylHinge(rho, l, r, sphereRad, hingePos, Nx, Ntheta);  
            
            testCase.verifyEqual(fb.Length, l);
            testCase.verifyEqual(fb.Radius, r);
            testCase.verifyEqual(fb.Cg, [0 0 0]);
            testCase.verifyEqual(fb.Zpos, 0);

            testCase.verifyEqual(size(fb.M), [8 8]);
            testCase.verifyEqual(size(fb.Dpto), [8 8]);
            testCase.verifyEqual(size(fb.K), [8 8]);

            sect = fb.WaterPlaneSec;

            figure;
            subplot(2,1,1);
            plot(sect(:,1), sect(:,2));
            title({'FloatingBodyUT - testFloatSphereEndCylHinge', ' ', 'Hinged Cylinder waterplane section'});
            axis equal;
            set(gca, 'xlim', [-l l], 'ylim', [-l l]);

            geo = fb.PanelGeo;

            subplot(2,1,2)
            plot(geo,'ShowSym');
            title('Hinged Cylinder panel geometry');
            axis equal;
        end
    end
end