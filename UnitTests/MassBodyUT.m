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
classdef MassBodyUT <  matlab.unittest.TestCase
    
    methods (Test)
        function testCyl(testCase)
            rho = 1;
            L = 10;
            R = 2;
            Nx = 40;
            Nr = 20;
            Ntheta = 20;
            cyl = makeMass_horCyl(rho, L, R, Nx, Nr, Ntheta);
            
            mex = pi*R^2*L;
            mac = cyl.Mass;
            
            testCase.verifyEqual(mac, mex, 'AbsTol', 1e-6);
            
            Mex = zeros(6,6);
            Mex(1,1) = mex;
            Mex(2,2) = mex;
            Mex(3,3) = mex;
            
            Mex(4,4) = mex*R^2/2;
            Mex(5,5) = 1/12*mex*(3*R^2 + L^2);
            Mex(6,6) = Mex(5,5);
            
            Mac = cyl.MassMatrix;
            
            testCase.verifyEqual(Mac, Mex, 'AbsTol', 1);
        end
        
        function testCylHinge(testCase)
            rho = 1;
            L = 10;
            R = 2;
            Nx = 40;
            Nr = 20;
            Ntheta = 20;
            cyl = makeMass_horCyl(rho, L, R, Nx, Nr, Ntheta);
            
            modes = ModesOfMotion();
            modes.Generalized = 1;
            % default hinge is about {x,z} = {0,0}
            modes.GenMotFuncs = HingeYFunc();
            
            cyl.Modes = modes;
            
            m = pi*R^2*L;
            
            Mex = zeros(6,6);
            Mex(1,1) = m;
            Mex(2,2) = m;
            Mex(3,3) = m;
            % TODO: add calculation to 'manual'
            Mex(4,4) = m*R^2/2;
            Mex(5,5) = 1/12*m*(3*R^2 + L^2);
            Mex(6,6) = Mex(5,5);
            Mex(3,7) = 1/4*m*L;
            Mex(7,3) = Mex(3,7);
            Mex(7,7) = Mex(5,5);
            
            Mac = cyl.MassMatrix;
            
            testCase.verifyEqual(Mac, Mex, 'AbsTol', 1);
        end
        
        function testCylHinge2(testCase)
            rho = 1;
            L = 12;
            R = 2;
            Nx = 60;
            Nr = 20;
            Ntheta = 20;
            cyl = makeMass_horCyl(rho, L, R, Nx, Nr, Ntheta);
            
            modes = ModesOfMotion();
            modes.Generalized = [1 1];
            hinge1 = HingeYFunc();
            xm = -2;
            hinge1.HingePos = [xm 0];
  
            hinge2 = HingeYFunc();
            xn = 2;
            hinge2.HingePos = [xn 0];
            modes.GenMotFuncs = [hinge1 hinge2];
            
            cyl.Modes = modes;
            
            m = pi*R^2*L;
            
            m = pi*R^2*L;

            Mex = zeros(8,8);
            Mex(1,1) = m;
            Mex(2,2) = m;
            Mex(3,3) = m;
            Mex(4,4) = 1/2*m*R^2;
            Mex(5,5) = 1/12*m*(3*R^2 + L^2);
            Mex(6,6) = Mex(5,5);

            Mex(7,7) = 1/12*m*(3*R^2 + L^2 + 12*xm^2);
            Mex(1,7) = -8/(3*pi)*m*R*(xm/L);
            Mex(7,1) = Mex(1,7);
            Mex(3,7) = 1/4*m*L*(1 + 4*(xm/L)^2);
            Mex(7,3) = Mex(3,7);
            Mex(5,7) = 1/12*m*(xm/L)*(4*xm^2 - 3*L^2);
            Mex(7,5) = Mex(5,7);

            Mex(8,8) = 1/12*m*(3*R^2 + L^2 + 12*xn^2);
            Mex(1,8) = -8/(3*pi)*m*R*(xn/L);
            Mex(8,1) = Mex(1,8);
            Mex(3,8) = 1/4*m*L*(1 + 4*(xn/L)^2);
            Mex(8,3) = Mex(3,8);
            Mex(5,8) =  1/12*m*(xn/L)*(4*xn^2 - 3*L^2);
            Mex(8,5) = Mex(5,8);
            Mex(7,8) = 1/12*m*(3*R^2 + L^2 + 12*xn*xm - 6*R^2*abs(xm - xn)/L + 4*abs(xm - xn)^3/L);     
            Mex(8,7) = Mex(7,8);
            
            Mac = cyl.MassMatrix;
            
            testCase.verifyEqual(Mac, Mex, 'AbsTol', 1);
        end
        
        function testSpheroid(testCase)
            rho = 1;
            L = 5;
            Ls = L/2;
            R = 1;
            Nx = 20;
            Nr = 20;
            Ntheta = 20;
            cyl = makeMass_spheroidEndHorCyl(rho, L, R, Ls, Nx, Nr, Ntheta);

            mex = 4/3*pi*R^2*Ls;

            mac = cyl.Mass;
            
            testCase.verifyEqual(mac, mex, 'RelTol', 1e-2);
        end
    end
end