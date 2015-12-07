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
classdef ConstraintMatCompUT < matlab.unittest.TestCase

    methods (Test)
        
        function test1(testCase)
            % one body, no hinge, global cg at cg 1
            cgs = [0, 0, 0];
            hin = [];
            
            P = ConstraintMatComp.HingedBodies(cgs, hin);
            
            Pex = diag([1 1 1 1 1 1]);
            
            for m = 1:6
                for n = 1:6
                    testCase.verifyEqual(P(m,n), Pex(m,n), 'AbsTol', 1e-12);
                end
            end
        end
        
        function test2(testCase)
            % one body, no hinge, global cg somewhere else
            cgs = [0, 0, 0];
            hin = [];
            org = [1, 0, 0];
            
            P = ConstraintMatComp.HingedBodies(cgs, hin, 'Origin', org);
            
            % composite body moves by x = 1, y = 2, z = 3
            % and rolls 1, pitches 2, and yaws 0.
            s = [1 2 3 1 2 0]';
            
            % new origin doesn't change wrt roll, yaw is 0, so only effect
            % is in pitch. Positive pitch of 2, inceases z position by 2*-(-1),
            % so z = 3 + 2*1, x and y are the same
            % anges of body 1 should all be the same
            qex = [1 2 5 1 2 0].';
            q = P.'*s;

            for m = 1:6
                testCase.verifyEqual(q(m), qex(m), 'AbsTol', 1e-12);
            end
            
            % composite body moves by x = 1, y = 2, z = 3
            % and rolls 1, pitches 0, and yaws 2.
            s = [1 2 3 1 0 2]';
            
            % new origin doesn't change wrt roll, pitch is 0, so only effect
            % is in yaw. Positive yaw of 2, decreases y position by 2*1,
            % so y = 2 - 2*1, x and z are the same
            % anges of body 1 should all be the same
            qex = [1 0 3 1 0 2].';
            q = P.'*s;

            for m = 1:6
                testCase.verifyEqual(q(m), qex(m), 'AbsTol', 1e-12);
            end
            
            % try another location
            hin = [];
            org = [1, -3, -2];
            
            P = ConstraintMatComp.HingedBodies(cgs, hin, 'Origin', org);
            
            % composite body moves by x = 1, y = 2, z = 3
            % and rolls 1, pitches 2, and yaws 3.
            s = [1 2 3 1 2 3]';
            
            % translation:  x = 1,      y = 2,      z = 3
            % roll:         x = +0,     y = -1*2,   z = +1*3
            % pitch:        x = +2*2,   y = +0,     z = +2*1
            % yaw:          x = -3*3,   y = -3*1,   z = +0
            % total:        x = -4,     y = -3,      z = 8
            % anges of body 1 should all be the same
            qex = [-4 -3 8 1 2 3].';
            q = P.'*s;

            for m = 1:6
                testCase.verifyEqual(q(m), qex(m), 'AbsTol', 1e-12);
            end
        end
        
        function test3(testCase)
            % two bodies, one hinge, 
            % global cg at body 1
            cgs = [-1, 0, -1; 1, 0, -2];
            hin = [0, 0, 0];
            
            P = ConstraintMatComp.HingedBodies(cgs, hin);

        end
        
        function test4(testCase)
            % two bodies, one hinge, 
            % global cg at hinge
            cgs = [-1, 0, 2; 3, 0, -4];
            hin = [0, 0, 0];
            org = hin;
            
            P = ConstraintMatComp.HingedBodies(cgs, hin, 'Origin', org);

            % composite body moves by x = 1, y = 2, z = 3
            % and rolls 1, pitches 2, yaws 3, and flexes 4
            s = [1 2 3 1 2 3 4]';
            
            qex = zeros(12, 1);
            % Body 1
            % translation:  x = 1,      y = 2,      z = 3
            % roll:         x = +0,     y = -1*2,   z = -1*0
            % pitch:        x = +2*2,   y = +0,     z = +2*1
            % yaw:          x = -3*0,   y = -3*1,   z = +0
            % total:        x = 5,      y = -3,      z = 5
            % anges of body 1 should all be the same
            qex(1) = 5;
            qex(2) = -3;
            qex(3) = 5;
            qex(4) = 1;
            qex(5) = 2;
            qex(6) = 3;
            
            % Body 2
            % pitch angle is: 4 + 2
            % translation:  x = 1,          y = 2,      z = 3
            % roll:         x = +0,         y = 1*4,    z = -1*0
            % pitch:        x = -(4+2)*4,   y = +0,     z = -(4+2)*3
            % yaw:          x = -3*0,       y = 3*3,   z = +0
            % total:        x = -23,        y = 15,     z = -15
            % anges of body 2 should all be the same
            qex(7) = -23;
            qex(8) = 15;
            qex(9) = -15;
            qex(10) = 1;
            qex(11) = 6;
            qex(12) = 3;
            
            q = P.'*s;

            for m = 1:12
                testCase.verifyEqual(q(m), qex(m), 'AbsTol', 1e-12);
            end
        end
        
        function test5(testCase)
            % two bodies, one hinge, 
            % global cg somewhere else
            P = ConstraintMatComp.Hinge(cgs, hin);
        end
        
         function test6(testCase)
            % three bodies, two hinges, 
            % global cg at cg 1
            P = ConstraintMatComp.Hinge(cgs, hin);
         end
        
         function test7(testCase)
            % three bodies, two hinges, 
            % global cg somewhere else
            P = ConstraintMatComp.Hinge(cgs, hin);
        end
       
        function test8(testCase)
            
            %             ______
            %    _______ |      |
            %   |   1   ||  2   |
            %   |_______||      |
            %            |______|
            %
            
            len1 = 5;
            wid1 = 4;
            hei1 = 2;
            
            len2 = len1;
            wid2 = wid1;
            hei2 = hei1;
            
            M1 = ConstraintMatCompUT.massBlock(len1, wid1, hei1);
            M2 = ConstraintMatCompUT.massBlock(len2, wid2, hei2);
            
            Mq = [M1; M2];
                        
            cg1 = [-len1/2, 0, 0];
            cg2 = [len2/2, 0, 0];
            
            cgs = [cg1; cg2];
            hin = [0 0 0];
            
            P = ConstraintMatComp.Hinge(cgs, hin);
            
            Ms = P*Mq*P.';
            
            Mex66 = ConstraintMatCompUT.massBlock(len1 + len2, wid1 + wid2, hei1 + hei2);
            
            Mex = zeros(7, 7);
            Mex(1:6, 1:6) = Mex66;
            Mex(7,7) = Mex66(5,5);
            
            for m = 1:7
                for n = 1:7
                    testCase.verifyEqual(Ms(m,n), Mex(m,n), 'AbsTol', 1e-12);
                end
            end
        end
    end
    
    methods (Static, Access = private)
        
        function [M] = massBlock(len, wid, hei)
            
            m = len*wid*hei;
            Ixx = m/12*(wid^2 + hei^2);
            Iyy = m/12*(len^2 + hei^2);
            Izz = m/12*(len^2 + wid^2);
            
            M = zeros(6, 6);
            M(1,1) = m;
            M(2,2) = m;
            M(3,3) = m;
            M(4,4) = Ixx;
            M(5,5) = Iyy;
            M(6,6) = Izz;
        end
    end
end