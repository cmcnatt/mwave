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
classdef HydroBodyUT < matlab.unittest.TestCase
    % Unit test for HydroBody class

    methods (Test)
        function testComputation(testCase)
            % Create the cylinder hydrobody

            load([mwavePath '\UnitTests\files\cyl6dof_WamOutUT']);

            cylHydBod = computeHydroBody(waveC, hydroF, floatB);

            load([mwavePath '\UnitTests\files\cyl1_verUT']);

            wComp = wComp1cyl;
            waveW = wWave1cyl;

            waveW.BodyMotions = wComp.Motions;

            pwr10 = -1;

            hComp = HydroBodyComp(cylHydBod,  wComp.IncWaves);            
            
            A = round2val(squeeze(hComp.A), pwr10);
            B = round2val(squeeze(hComp.B), pwr10);
            Fex = round2val(squeeze(hComp.Fex).', pwr10);

            Aw = round2val(squeeze(wComp.A), pwr10);
            Bw = round2val(squeeze(wComp.B), pwr10);
            Fexw = round2val(squeeze(wComp.Fex).', pwr10);
            
%             errA = matErr(A, Aw);
%             errB = matErr(B, Bw);
%             errFex = matErr(Fex, Fexw);
            
            testCase.verifyEqual(size(A), [6 6]);
            testCase.verifyEqual(size(B), [6 6]);
            testCase.verifyEqual(size(Fex), [6 5]);
            
            for m = 1:6
                testCase.verifyEqual(Fex(m), Fexw(m), 'AbsTol', 10^pwr10);
                for n = 1:6
                    testCase.verifyEqual(A(m,n), Aw(m,n), 'AbsTol', 10^pwr10);
                    testCase.verifyEqual(B(m,n), Bw(m,n), 'AbsTol', 10^pwr10);
                end
            end
        end
    end
    
end