classdef ViscDampCoefUT < matlab.unittest.TestCase
    
    methods (Test)
        function areaBox(testCase)
            
            x = [-1 1];
            y = [-2 2];
            z = [-3 3];
            centRot = [0 0 0];
            
            A = ViscDampCoef.ComputeAreaFromBox(x, y, z, centRot);
            
            Aexp = [4*6, 2*6, 2*4, ...
                2/4*2*(2^2+3^2)^2, ...
                4/4*2*(1^2+3^2)^2, ...
                6/4*2*(1^2+2^2)^2];
            
            testCase.assertEqual(A, Aexp, 'relTol', 10^-10);
            
            centRot = [-1 2 -0.5];
            x = x + centRot(1);
            y = y + centRot(2);
            z = z + centRot(3);
            
            A = ViscDampCoef.ComputeAreaFromBox(x, y, z, centRot);
            
            testCase.assertEqual(A, Aexp, 'relTol', 10^-10);
        end
        
        function coefArea(testCase)
            
            x = [-1 1];
            y = [-2 2];
            z = [-3 3];
            centRot = [0 0 0];
            
            coefs = ViscDampCoef.ComputeAreaFromBox(x, y, z, centRot, 'makeCoefs');
            
            Aexp = [4*6, 2*6, 2*4, ...
                2/4*2*(2^2+3^2)^2, ...
                4/4*2*(1^2+3^2)^2, ...
                6/4*2*(1^2+2^2)^2];
            
            for n = 1:6
                testCase.assertEqual(coefs(n).A, Aexp(n), 'relTol', 10^-10);
            end
        end
    end
end