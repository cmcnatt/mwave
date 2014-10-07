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
classdef ModesOfMotionUT < matlab.unittest.TestCase

    methods (Test)
       
        function test1(testCase)
            modes = ModesOfMotion();

            testCase.verifyTrue(modes.Surge);
            
            testCase.verifyTrue(modes.Sway);
            
            testCase.verifyTrue(modes.Heave);
            
            testCase.verifyTrue(modes.Roll);
            
            testCase.verifyTrue(modes.Pitch);
            
            testCase.verifyTrue(modes.Yaw);

            testCase.verifyTrue(~modes.HasGen);
            
            testCase.verifyEqual(modes.DoF, 6);
        end
        
        function test2(testCase)
            vect = [1 0 1 1 0 1];
            modes = ModesOfMotion(vect);
            
            testCase.verifyTrue(modes.Surge);
            
            testCase.verifyTrue(~modes.Sway);
            
            testCase.verifyTrue(modes.Heave);
            
            testCase.verifyTrue(modes.Roll);
            
            testCase.verifyTrue(~modes.Pitch);
            
            testCase.verifyTrue(modes.Yaw);

            testCase.verifyTrue(~modes.HasGen);
            
            testCase.verifyEqual(modes.DoF, 4);
        end
        
        function test3(testCase)
            vect = [1 0 1 1 0 1 5];
            modes = ModesOfMotion(vect);
            
            testCase.verifyTrue(modes.Surge);
            
            testCase.verifyTrue(~modes.Sway);
            
            testCase.verifyTrue(modes.Heave);
            
            testCase.verifyTrue(modes.Roll);
            
            testCase.verifyTrue(~modes.Pitch);
            
            testCase.verifyTrue(modes.Yaw);
            
            testCase.verifyTrue(modes.HasGen);
            
            testCase.verifyEqual(modes.Generalized, 5);
            
            testCase.verifyEqual(modes.DoF, 9);
            
            motions = modes.Motions;
            
            testCase.verifyEqual(motions{1}, 'Surge');
            testCase.verifyEqual(motions{2}, 'Heave');
            testCase.verifyEqual(motions{3}, 'Roll');
            testCase.verifyEqual(motions{4}, 'Yaw');
            testCase.verifyEqual(motions{5}, 'Gen1');
            testCase.verifyEqual(motions{6}, 'Gen2');
            testCase.verifyEqual(motions{7}, 'Gen3');
            testCase.verifyEqual(motions{8}, 'Gen4');
            testCase.verifyEqual(motions{9}, 'Gen5');
        end
    end
end