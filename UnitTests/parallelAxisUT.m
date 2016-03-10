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
classdef parallelAxisUT < matlab.unittest.TestCase

    methods (Test)
       
        function test1(testCase)
            I = diag([1, 2, 3]);
            r = [-1, 2, -3];
            
            m = 5;
            
            Iout = parallelAxis(I, m, r);
            
            Iex = [1 + 5*(2^2 + 3^2), -5*(-1)*2, -5*(-1)*(-3);
                -5*(-1)*2, 2 + 5*(1^2 + 3^2), -5*2*(-3);
                -5*(-1)*(-3), -5*2*(-3), 3 + 5*(1^2 + 2^2)];
            
            for m = 1:3
                for n = 1:3
                    testCase.verifyEqual(Iout(m,n), Iex(m,n), ...
                        'AbsTol', 1e-12);
                end
            end
        end
    end
end