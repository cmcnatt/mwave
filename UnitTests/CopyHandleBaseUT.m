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
classdef CopyHandleBaseUT < matlab.unittest.TestCase

    methods (Test)
        function test1(testCase)
            base = CopyHandleTest(2.1, 2.2);
            base.PubProp1 = 1;
            base.PubProp2 = 2;
            
            base.PriProp1 = 1.1;
            base.ProProp1 = 1.2;
            
            deep = copy(base);
            deep2 = base.copy;
            
            shallow = base;
            
            % change base values
            base.PubProp1 = 10;
            base.PubProp2 = 20;
            base.PriProp1 = 11;
            base.ProProp1 = 12;
            base.SetPriProp2(21);
            base.SetProProp2(22);
            
            testCase.assertEqual(shallow.PubProp1, base.PubProp1);
            testCase.assertEqual(shallow.PubProp2, base.PubProp2);
            testCase.assertEqual(shallow.PriProp1, base.PriProp1);
            testCase.assertEqual(shallow.GetPriProp2, base.GetPriProp2);
            testCase.assertEqual(shallow.ProProp1, base.ProProp1);
            testCase.assertEqual(shallow.GetProProp2, base.GetProProp2);
            
            testCase.assertEqual(deep.PubProp1, 1);
            testCase.assertEqual(deep.PubProp2, 2);
            testCase.assertEqual(deep.PriProp1, 1.1);
            testCase.assertEqual(deep.GetPriProp2, 2.1);
            testCase.assertEqual(deep.ProProp1, 1.2);
            testCase.assertEqual(deep.GetProProp2, 2.2);
            
            testCase.assertEqual(deep2.PubProp1, 1);
            testCase.assertEqual(deep2.PubProp2, 2);
            testCase.assertEqual(deep2.PriProp1, 1.1);
            testCase.assertEqual(deep2.GetPriProp2, 2.1);
            testCase.assertEqual(deep2.ProProp1, 1.2);
            testCase.assertEqual(deep2.GetProProp2, 2.2);
        end
    end
end