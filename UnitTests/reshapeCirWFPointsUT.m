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
classdef reshapeCirWFPointsUT < matlab.unittest.TestCase

    methods (Test)
       
        function test1(testCase)
            r = 1;
            Ntheta = 4;
            zin = [0 -2 -4]';

            points = makeCirWFPoints(r, Ntheta, zin);

            N = size(points,1);
            eta0 = 1:1:N;

            %  r  |  z  |  theta  | eta
            %----------------------------
            %  1  |  0  |    0    |  1
            %  1  |  0  |   pi/2  |  2
            %  1  |  0  |   pi    |  3
            %  1  |  0  |  3pi/2  |  4
            %  1  | -2  |    0    |  5
            %  1  | -2  |   pi/2  |  6
            %  1  | -2  |   pi    |  7
            %  1  | -2  |  3pi/2  |  8
            %  1  | -4  |    0    |  9
            %  1  | -4  |   pi/2  |  10
            %  1  | -4  |   pi    |  11
            %  1  | -4  |  3pi/2  |  12
            
            [r, theta, z, eta] = reshapeCirWFPoints(points, eta0);
            
            testCase.verifyEqual(r, 1);
            testCase.verifyEqual(theta, [0, pi/2, pi, 3*pi/2]);
            testCase.verifyEqual(z, zin);

            fail = false;

            i = 0;
            for m = 1:2
                for n = 1:4;
                    i = i+1;
                    if (eta(m,n) ~= i)
                        fail = true;
                        break;
                    end
                end
            end
            
            testCase.verifyTrue(~fail);
        end
    end
end