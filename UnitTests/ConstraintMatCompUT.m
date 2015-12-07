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