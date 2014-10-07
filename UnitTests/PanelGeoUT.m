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
classdef PanelGeoUT < matlab.unittest.TestCase

    methods (Test)
       
        function testRotation(testCase)
            
            len = 10;
            rad = 0.5;
            
            Nx = 120;
            Ntheta = 16;
            
            sphereRad = 0.1*len;
            hingePos = [-len/6; len/6];
            
            panelGeo = makePanel_spheroidEndHorCylHinge(len, rad, sphereRad, hingePos, Nx, Ntheta);
            
            figure;
            subplot(3,1,1);
            plot(panelGeo);
            axis equal;
            view([0 0]);
            title({'PanelGeoUT - testRotation', 'Unrotated body'});
            
            % pitch
            ax = [0 1 0];
            angle = 30/180*pi;
            
            panelGeo.Rotate(ax, angle);
            
            subplot(3,1,2);
            plot(panelGeo);
            axis equal;
            view([0 0]);
            title({'Pitched 30 degrees'});
            
            % rotate about hinge
            ax = [0 1 0];
            angle = 30/180*pi;
            org = [hingePos(1), 0, 0];
            
            hingeGeo = makePanel_spheroidEndHorCylHinge(len, rad, sphereRad, hingePos, Nx, Ntheta);
            
            hingeGeo.Rotate(ax, angle, 'Origin', org);
            
            subplot(3,1,3);
            plot(hingeGeo);
            axis equal;
            view([0 0]);
            title({'Rot 30 deg about front hinge'});
        end
        
        function testSplit(testCase)
            
            len = 10;
            rad = 0.5;
            
            Nx = 120;
            Ntheta = 16;
            
            sphereRad = 0.1*len;
            hingePos = [-len/6; len/6];
            
            panelGeo = makePanel_spheroidEndHorCylHinge(len, rad, sphereRad, hingePos, Nx, Ntheta);
            
            figure;
            subplot(3,1,1);
            plot(panelGeo);
            axis equal;
            %view([0 0]);
            title({'PanelGeoUT - testSplit', 'Whole body'});
            
            % split
            
            planePoint = [hingePos(1), 0, 0];
            planeNorm = [1 0 0];
            
            [backGeo, frontGeo] = panelGeo.Split(planePoint, planeNorm);
            
            backGeo.Translate([-1 0 0]);
            frontGeo.Translate([1 0 0]);
            
            subplot(3,1,2);
            plot(backGeo);
            plot(frontGeo);
            axis equal;
            %view([0 0]);
            title({'Split at front hinge'});
            
            % split then rotate
            backGeo.Translate([1 0 0]);
            frontGeo.Translate([-1 0 0]);

            ax = [0 1 0];
            angle = 15/180*pi;
            org = [hingePos(1), 0, 0];
            
            backGeo.Rotate(ax, angle, 'Origin', org);
            frontGeo.Rotate(ax, -angle, 'Origin', org);
            
            
            subplot(3,1,3);
            surf(backGeo);
            surf(frontGeo);
            axis equal;
            title({'Split then rotated 30 deg about front hinge'});
            
            
        end
    end
end