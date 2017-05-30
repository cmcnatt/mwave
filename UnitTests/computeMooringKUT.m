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
classdef computeMooringKUT < matlab.unittest.TestCase

    methods (Test)
        
        function test0(testCase)
            % nlin = 1
            knum = 8;
            ks = knum/(sqrt(2)/2);
            cg = [0 0 0];
            pos = [2 0 -1];
            anghs = 0;
            angvs = pi/4;
           
            K = computeMooringK(cg, pos, ks, anghs, angvs);
            
            Kexp = zeros(6,6);
            Kexp(1,1) = knum;
            Kexp(3,3) = knum;
            
            r = pos - cg;
            
            Kexp(3,5) = r(1)*knum;
            Kexp(5,3) = r(1)*knum;
            Kexp(5,5) = r(1)^2*knum + r(3)^2*knum;
            Kexp(1,5) = -r(3)*knum;Kexp(1,5) = -r(3)*knum;
            Kexp(5,1) = -r(3)*knum;
            
            for m = 1:6
                for n = 1:6
                    testCase.verifyEqual(K(m,n), Kexp(m,n), 'AbsTol', 1e-9);
                end
            end
            
            x = [1 0 1 0 1 0].';
            
            Fk = -K*x;
            
            Tk = -r(3)*Fk(1) + r(1)*Fk(3);
            
            testCase.verifyEqual(Tk, Fk(5), 'AbsTol', 1e-11);
        end
       
        function test1(testCase)
            % pos at cg
            % single vertical mooring line
            % no surge, sway, roll, pitch, yaw
            % just heave
            
            % nlin = 1
            ks = 8;
            cg = [1 1 1];
            pos = cg;
            anghs = 0;
            angvs = pi/2;
           
            K = computeMooringK(cg, pos, ks, anghs, angvs);
            
            Kexp = zeros(6,6);
            Kexp(3,3) = ks;
            
            for n = 1:6
                testCase.verifyEqual(K(n,n), Kexp(n,n), 'AbsTol', 1e-9);
            end
        end
        
        function test2(testCase)
            % pos at cg
            % two vertical mooring lines going forward and backward
            % at 45 deg (pi/4)
            % no sway, roll, pitch, yaw
            % just surge, heave
            
            % nlin = 2
            ks = [8 4];
            cg = [1 1 1];
            pos = cg;
            anghs = [0 -pi];
            angvs = [pi/4 pi/4];
           
            K = computeMooringK(cg, pos, ks, anghs, angvs);
            
            Kexp = zeros(6,6);
            Kexp(1,1) = sqrt(2)/2*(ks(1)+ks(2));
            Kexp(3,3) = Kexp(1,1);
            
            for n = 1:6
                testCase.verifyEqual(K(n,n), Kexp(n,n), 'AbsTol', 1e-9);
            end
        end
        
        function test3(testCase)
            % pos at cg
            % two vertical mooring lines going side to side
            % at 45 deg (pi/4) from horizontal
            % no surge, roll, pitch, yaw
            % just sway, heave
            
            % nlin = 2
            ks = [8 4];
            cg = [1 1 1];
            pos = cg;
            anghs = [-pi/2 pi/2];
            angvs = [pi/4 pi/4];
           
            K = computeMooringK(cg, pos, ks, anghs, angvs);
            
            Kexp = zeros(6,6);
            Kexp(2,2) = sqrt(2)/2*(ks(1)+ks(2));
            Kexp(3,3) = Kexp(2,2);
            
            for n = 1:6
                testCase.verifyEqual(K(n,n), Kexp(n,n), 'AbsTol', 1e-9);
            end
        end
        
        function test4(testCase)
            % zpos = zcg, r = 2
            % one vertical mooring line
            % no surge, sway, roll, yaw
            % just heave, pitch
            
            % nlin = 1
            ks = 8;
            cg = [1 1 1];
            r = 2;
            pos = [cg(1) + r, cg(2), cg(3)];
            anghs = 0;
            angvs = pi/2;
           
            K = computeMooringK(cg, pos, ks, anghs, angvs);
            
            Kexp = zeros(6,6);
            Kexp(3,3) = ks;
            Kexp(5,5) = ks*r^2;
            
            for n = 1:6
                testCase.verifyEqual(K(n,n), Kexp(n,n), 'AbsTol', 1e-9);
            end
        end
        
        function test5(testCase)
            % zpos = zcg, r = 2
            % two vertical mooring lines at two positions
            % no surge, sway, roll, yaw
            % just heave, pitch
            
            % nlin = 1
            ks = [8, 4];
            cg = [1 1 1];
            r = 2;
            pos = [cg(1) - r, cg(2), cg(3); cg(1) + r, cg(2), cg(3)];
            anghs = [0 0];
            angvs = [pi/2 pi/2];
           
            K = computeMooringK(cg, pos, ks, anghs, angvs);
            
            Kexp = zeros(6,6);
            Kexp(3,3) = sum(ks);
            Kexp(5,5) = sum(ks)*r^2;
            
            for n = 1:6
                testCase.verifyEqual(K(n,n), Kexp(n,n), 'AbsTol', 1e-9);
            end
        end
        
        function test6(testCase)
            % zpos = zcg, r = 2
            % two vertical mooring lines at two positions
            % no surge, sway, ptich, yaw
            % just heave, roll
            
            % nlin = 1
            ks = [8, 4];
            cg = [1 1 1];
            r = 2;
            pos = [cg(1), cg(2)-r, cg(3); cg(1), cg(2)+r, cg(3)];
            anghs = [0 0];
            angvs = [pi/2 pi/2];
           
            K = computeMooringK(cg, pos, ks, anghs, angvs);
            
            Kexp = zeros(6,6);
            Kexp(3,3) = sum(ks);
            Kexp(4,4) = sum(ks)*r^2;
            
            for n = 1:6
                testCase.verifyEqual(K(n,n), Kexp(n,n), 'AbsTol', 1e-9);
            end
        end
        
        function test7(testCase)
            % rx = 3, rz = 2 
            % two mooring lines at one position
            % at 45 deg
            % no sway, no roll, no yaw
            
            % nlin = 1
            ks = [8, 4];
            cg = [1 1 1];
            rx = 3;
            rz = 2;
            pos = [cg(1)+rx, cg(2), cg(3)-rz];
            anghs = [0 pi];
            angvs = [pi/4 pi/4];
           
            K = computeMooringK(cg, pos, ks, anghs, angvs);
            
            Kexp = zeros(6,6);
            Kexp(1,1) = sqrt(2)/2*sum(ks);
            Kexp(3,3) = sqrt(2)/2*sum(ks);
            Kexp(5,5) = sqrt(2)/2*sum(ks)*(rx^2 + rz^2);
            
            for n = 1:6
                testCase.verifyEqual(K(n,n), Kexp(n,n), 'AbsTol', 1e-9);
            end
        end
        
        function test8(testCase)
            % ry = 3, rz = 2 
            % two mooring lines at one position
            % at 30 deg
            % no sway, no roll, no yaw
            
            % nlin = 1
            ks = [8, 4];
            cg = [1 1 1];
            ry = 3;
            rz = 2;
            pos = [cg(1), cg(2)+ry, cg(3)-rz];
            anghs = [-pi/2 pi/2];
            angvs = [pi/6 pi/6];
           
            K = computeMooringK(cg, pos, ks, anghs, angvs);
            
            Kexp = zeros(6,6);
            Kexp(2,2) = sqrt(3)/2*sum(ks);
            Kexp(3,3) = 1/2*sum(ks);
            Kexp(4,4) = (sqrt(3)/2*rz^2 + 1/2*ry^2)*sum(ks);
            
            for n = 1:6
                testCase.verifyEqual(K(n,n), Kexp(n,n), 'AbsTol', 1e-9);
            end
        end
        
        function test9(testCase)
            % ry = 3, rz = 2 
            % two mooring lines at one position
            % at 30 deg
            % no sway, no roll, no yaw
            
            % nlin = 1
            ks = [8, 4];
            cg = [1 1 1];
            ry = 3;
            rz = 2;
            pos = [cg(1), cg(2)+ry, cg(3)-rz];
            anghs = [-pi/2 pi/2];
            angvs = [pi/6 pi/6];
           
            K = computeMooringK(cg, pos, ks, anghs, angvs);
            
            Kexp = zeros(6,6);
            Kexp(2,2) = sqrt(3)/2*sum(ks);
            Kexp(3,3) = 1/2*sum(ks);
            Kexp(4,4) = (sqrt(3)/2*rz^2 + 1/2*ry^2)*sum(ks);
            
            for n = 1:6
                testCase.verifyEqual(K(n,n), Kexp(n,n), 'AbsTol', 1e-9);
            end
        end
        
        function test10(testCase)
            % three lines: 0, 2*pi/3, -2*pi/3
            % hor: pi/6
            % k = 8
            
            % nlin = 1
            ks = [8, 8, 8];
            cg = [0 0 0];
            pos = [6, 0, -1];
            anghs = [0, -2*pi/3, 2*pi/3];
            angvs = [pi/6 pi/6 pi/6];
           
            K = computeMooringK(cg, pos, ks, anghs, angvs);
            
            Kexp = zeros(6,6);
            Kexp(1,1) = sqrt(3)/2*(1*ks(1) + 1/2*ks(2) + 1/2*ks(3));
            Kexp(2,2) = sqrt(3)/2*(sqrt(3)/2*ks(2) + sqrt(3)/2*ks(3));
            Kexp(3,3) = 1/2*sum(ks);
            Kexp(4,4) = Kexp(2,2)*abs(pos(3))^2;
            Kexp(5,5) = Kexp(1,1)*abs(pos(3))^2 + Kexp(3,3)*pos(1)^2;
            Kexp(6,6) = Kexp(2,2)*pos(1)^2;
            
            for n = 1:6
                testCase.verifyEqual(K(n,n), Kexp(n,n), 'AbsTol', 1e-9);
            end
        end
    end
end