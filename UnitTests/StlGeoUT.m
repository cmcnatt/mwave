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
            
            filePath = [mwavePath '\UnitTests\files\'];
            fwd = StlGeo;
            fwd.Read([filePath 'forward.stl']);
            
            aft = StlGeo;
            aft.Read([filePath 'aft.stl']);
            
            figure;
            fwd.Plot;
            %aft.Plot;
            axis equal
        end
    end
end