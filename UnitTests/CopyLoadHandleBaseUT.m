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
classdef CopyLoadHandleBaseUT < matlab.unittest.TestCase

    methods (Test)
        function testCopy1(testCase)
            base = CopyLoadHandleTest(2.1, 2.2);
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
        
        function testLoad1(testCase)
            
            pub1 = 1;
            pub2 = 2;
            pri11 = 1.1;
            pri12 = 1.2;
            pri21 = 2.1;
            pri22 = 2.2;
            
            
            base = CopyLoadHandleTest(pri21, pri22);
            base.PubProp1 = pub1;
            base.PubProp2 = pub2;
            
            base.PriProp1 = pri11;
            base.ProProp1 = pri12;
            
            save([matTempPath '\loadHandleTest1'], 'base');
            
            clear base;
            
            load([matTempPath '\loadHandleTest1'])
            
            testCase.assertEqual(base.PubProp1, pub1);
            testCase.assertEqual(base.PubProp2, pub2);
            testCase.assertEqual(base.PriProp1, pri11);
            testCase.assertEqual(base.GetPriProp2, pri21);
            testCase.assertEqual(base.ProProp1, pri12);
            testCase.assertEqual(base.GetProProp2, pri22);
        end
        
        function testLoad2(testCase)
            pub1 = 1;
            pub2 = 2;
            pri11 = 1.1;
            pri12 = 1.2;
            pri21 = 2.1;
            pri22 = 2.2;
            
            pub3 = 3;
            pri13 = 1.3;
            pri23 = 2.3;
            
            load([myMatPath '\mwave\UnitTests\files\loadHandleTest2'])
            
            testCase.assertEqual(base2.PubProp1, pub1);
            testCase.assertEqual(base2.PubProp2, pub2);
            testCase.assertEqual(base2.PriProp1, pri11);
            testCase.assertEqual(base2.GetPriProp2, pri21);
            testCase.assertEqual(base2.ProProp1, pri12);
            testCase.assertEqual(base2.GetProProp2, pri22);
            
            unSetProps = base2.GetUnsetProps;
            
            testCase.assertEqual(unSetProps.PubProp3, pub3);
            testCase.assertEqual(unSetProps.priProp3, pri13);
            testCase.assertEqual(unSetProps.proProp3, pri23);
        end
    end
    
    methods (Static)
        
        function [] = SaveExampleMatFile()
            pub1 = 1;
            pub2 = 2;
            pri11 = 1.1;
            pri12 = 1.2;
            pri21 = 2.1;
            pri22 = 2.2;
            
            pub3 = 3;
            pri13 = 1.3;
            pri23 = 2.3;

            % create a version of CopyHandleTest that has these properties,
            % save it, and then delete the properties

            base2 = CopyLoadHandleTest(pri21, pri22);
            base2.PubProp1 = pub1;
            base2.PubProp2 = pub2;
            base2.PubProp3 = pub3;

            base2.PriProp1 = pri11;
            base2.ProProp1 = pri12;
            base2.PriProp3 = pri13;
            base2.ProProp3 = pri23;

            save([myMatPath '\mwave\UnitTests\files\loadHandleTest2'], 'base2');
        end
    end
end