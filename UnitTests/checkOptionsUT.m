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
classdef checkOptionsUT < matlab.unittest.TestCase

    methods (Test)
       
        function test1(testCase)
            
            varargin = {'ShowMovie', 'Why', 'k', 5, 0.1, 'Radiated', 'Mlim', 20};
            
            posOpts = {{'ShowMovie'}, {'Mlim', 1}, {'Radiated'}, {'Scattered'}, {'Why', 3}}; 
            
            [opts, args] = checkOptions(posOpts, varargin);
            
            showMov = opts(1);
            useMlim = opts(2);
            Mlim = args{2};
            rad = opts(3);
            scat = opts(4);
            why = opts(5);
            whyArgs = args{5};
            
            
            
            testCase.verifyTrue(showMov == true);
            testCase.verifyTrue(useMlim == true);
            testCase.verifyTrue(rad == true);
            testCase.verifyTrue(scat == false);
            testCase.verifyTrue(why == true);
            
            testCase.verifyTrue(Mlim == 20);
            testCase.verifyTrue(strcmp(whyArgs{1}, 'k'));
            testCase.verifyTrue(whyArgs{2} == 5);
            testCase.verifyTrue(whyArgs{3} == 0.1);
            
        end
        
    end
end