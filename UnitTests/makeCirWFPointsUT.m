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
classdef makeCirWFPointsUT < matlab.unittest.TestCase

    methods (Test)
       
        function test1(testCase)
            r = 1;
            Ntheta = 4;
            zin = [0 -2 -4]';

            points = makeCirWFPoints(r, Ntheta, zin);
            
            Np = size(points, 1);
            testCase.verifyEqual(Np, Ntheta*length(zin));
            
            rp = sqrt(points(:,1).^2 + points(:,2).^2);
            xy = [1 0; 0 1; -1 0; 0 -1];
            
            for n = 1:Np
                testCase.verifyEqual(rp(n), r, 'AbsTol', 1e-12);
                testCase.verifyEqual(points(n, 1), xy(round(mod(n, 4.1)), 1), 'AbsTol', 1e-10);
                testCase.verifyEqual(points(n, 2), xy(round(mod(n, 4.1)), 2), 'AbsTol', 1e-10);
                testCase.verifyEqual(points(n, 3), zin(ceil(n/4)), 'AbsTol', 1e-10);
            end
        end
        
    end
end