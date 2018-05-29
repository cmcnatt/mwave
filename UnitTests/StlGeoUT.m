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
classdef StlGeoUT <  matlab.unittest.TestCase
    
    methods (Test)
        function test1(testCase)            
            [fwd, aft] = StlGeoUT.getStl6DOFs();
            
            figure;
            plot(fwd);
            plot(aft);
            set(gca, 'view', [-16 3]);
            axis equal
        end
        
        function test2(testCase)            
            [fwd, aft] = StlGeoUT.getStl6DOFs();
            
            fwd.Position = [10 10 10];
            aft.Position = [10 10 10];
            
            fwd.Orientation = [0 5/180*pi 0];
            aft.Orientation = [0 -5/180*pi 0];
            
            figure;
            plot(fwd);
            plot(aft);
            set(gca, 'view', [-16 3]);
            axis equal
        end
        
        function test3(testCase)  
            
            pit = 15/180*pi;
            flex = -30/180*pi;
            [fwd, aft] = StlGeoUT.getStl6DOFs();
            
            fwd.Position = [0 0 0];
            aft.Position = [0 0 0];
            
            fwd.Orientation = [0 pit 0];
            aft.Orientation = [0 flex+pit 0];
            
            geo = StlGeoUT.getStlHinge();
            geo.Position = fwd.Position;
            geo.Orientation = fwd.Orientation;
            geo.FlexAngle = flex;
            
            figure;
            subplot(2,1,1);
            plot(fwd);
            plot(aft);
            box on
            set(gca, 'view', [-16 3]);
            
            axis equal
            
            subplot(2,1,2);
            plot(geo);
            axis equal;
            box on
            set(gca, 'view', [-16 3]);
        end
        
        function test4(testCase)  
            
            sur = 20;
            pit = 15/180*pi;
            flex = -30/180*pi;
            yaw = 45/180*pi;
            [fwd, aft] = StlGeoUT.getStl6DOFs();
            
            fwd.Position = [sur 0 0];
            aft.Position = [sur 0 0];
            
            fwd.Orientation = [0 pit yaw];
            aft.Orientation = [0 flex+pit yaw];
            
            geo = StlGeoUT.getStlHinge();
            geo.Position = fwd.Position;
            geo.Orientation = fwd.Orientation;
            geo.FlexAngle = flex;
            
            figure;
            subplot(2,1,1);
            plot(fwd);
            plot(aft);
            set(gca, 'view', [-16 3]);
            
            axis equal
            
            subplot(2,1,2);
            plot(geo);
            axis equal;
            set(gca, 'view', [-16 3]);
        end
        
        function test5(testCase)
            a = 2;
            omega = 2*pi/6;
            beta = 0;
            h = 50;
            geo = StlWave(a, omega, beta, h);
            
            geo.X = -100:5:100;
            geo.Y = -100:5:100;
            geo.Time = 0;
            
            figure;
            plot(geo);
            set(gca, 'view', [-16 3]);
            axis equal
            
        end
    end
    
    methods (Static)
        function [fwd, aft] = getStl6DOFs()
            
            filePath = [mwavePath '\UnitTests\files\'];
            fwd = Stl6DOFGeo;
            fwd.Read([filePath 'forward.stl']);
            
            aft = Stl6DOFGeo;
            aft.Read([filePath 'aft.stl']);
            
            fwd.Cg = 50*[-0.2727, 0, 0];             
            aft.Cg = 50*[0.5764, 0, -0.02];
        end
        
        function [geo] = getStlHinge()
            
            filePath = [mwavePath '\UnitTests\files\'];
            geo = StlHingeGeo;
            geo.ReadForward([filePath 'forward.stl']);
            geo.ReadAft([filePath 'aft.stl']);
            
            geo.ForwardCg = 50*[-0.2727, 0, 0];             
            geo.AftCg = 50*[0.5764, 0, -0.02];
            geo.HingeLoc = [0 0 0];
        end
    end
end