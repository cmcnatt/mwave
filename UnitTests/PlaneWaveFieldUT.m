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
classdef PlaneWaveFieldUT < matlab.unittest.TestCase

    methods (Test)
       
        function testElVel(testCase)
            rho = 1000;

            x = -10:0.5:10;
            y = -10:0.5:10;

            [X, Y] = meshgrid(x, y);

            wc = PlaneWaves(1, 2, 0, 10);

            wf = PlaneWaveField(rho, wc, 1, X, Y);

            eta = wf.Elevation;

            figure;
            subplot(2,1,1);
            pcolor(X,Y,real(eta{1}));
            shading flat;
            axis equal;
            axis tight;
            title({'PlaneWaveFieldUT - testElVel', ' ', 'Elevation'})
            
            vel = wf.Velocity;

            subplot(2,1,2);
            quiver(X,Y,real(vel{1}), real(vel{2}));
            axis equal;
            axis tight;
            title('Velocity');
        end
        
        function testNoVel(testCase)
            rho = 1000;

            x = -10:0.5:10;
            y = -10:0.5:10;

            [X, Y] = meshgrid(x, y);

            wc = PlaneWaves(1, 2, 0, 10);

            % default with velocity
            wf = PlaneWaveField(rho, wc, 1, X, Y);
            
            testCase.verifyTrue(wf.HasVelocity);
            
            try
                vel = wf.Velocity;
                % Should not cause an error
            catch
                testCase.verifyFail();
            end
            
            % default 2 with velocity#
            pnts = [0 0 0; 1 1 0];
            wf = PlaneWaveField(rho, wc, false, pnts);
            
            testCase.verifyTrue(wf.HasVelocity);
            
            try
                vel = wf.Velocity;
                % Should not cause an error
            catch
                testCase.verifyFail();
            end
            
            % no velocity case
            wf = PlaneWaveField(rho, wc, 1, X, Y, 'NoVel');
            
            testCase.verifyTrue(~wf.HasVelocity);
            
            try
                vel = wf.Velocity;
                % Should cause an error
                testCase.verifyFail();
            catch
                % Pass - does throw error
            end
            
            % no velocity case
            wf = PlaneWaveField(rho, wc, 0, pnts, 'NoVel');
            
            testCase.verifyTrue(~wf.HasVelocity);
            
            try
                vel = wf.Velocity;
                % Should cause an error
                testCase.verifyFail();
            catch
                % Pass - does throw error
            end
        end
        
        function testMultiWaveComp(testCase)
            rho = 1000;

            x = -10:0.1:10;
            y = -10:0.1:10;

            [X, Y] = meshgrid(x, y);

            T = [4, 3, 2, 1];
            nT = length(T);
            
            wcs = PlaneWaves(ones(size(T)), T, zeros(nT, 1), 10);

            wf = PlaneWaveField(rho, wcs, 1, X, Y);

            eta = wf.Elevation;

            figure;
            for n = 1:nT
                subplot(2,2,n);
                pcolor(X,Y,real(eta{n}));
                shading flat;
                axis equal;
                axis tight;
                if (n == 1)
                    title({'PlaneWaveFieldUT - testMultiWaveComp', ' ', ['T = ' num2str(T(n))]})
                else
                    title(['T = ' num2str(T(n))]);
                end
            end

            hs = wf.SigWaveHeight;

            T = [3, 2];
            beta = [0 pi/2];

            wcs = PlaneWaves([1 1], T, beta, 10);
            
            try
                wf = PlaneWaveField(rho, wcs, 1, X, Y);
                % Should cause error due to different directions
                testCase.verifyFail();
            catch
                % Pass if no error
            end
            
            wcs = PlaneWaves([1 1], [1 1], beta, 10);
             
             try
                wf = PlaneWaveField(rho, wcs, 1, X, Y);
                % Should cause error due to same frequency
                testCase.verifyFail();
            catch
                % Pass if no error
            end
        end
        
    end
end