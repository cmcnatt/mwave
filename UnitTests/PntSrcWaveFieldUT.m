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
classdef PntSrcWaveFieldUT < matlab.unittest.TestCase

    methods (Test)
       
        function testHeaveCenter(testCase)
            rho = 1000;

            x = -10:0.1:10;
            y = -10:0.1:10;

            [X, Y] = meshgrid(x, y);

            wc = PlaneWaves(1, 1, 0, 10);
            loc = [0 0];

            wf = PntSrcWaveField('Heave', loc, rho, wc, 1, X, Y);

            eta = wf.Elevation;
            eta = eta{1};

            t = 0:0.01:1;
            figure;

            subplot(2,1,1);
            for j = 1:length(t)
                cla(gca);
                pcolor(X,Y,real(eta*exp(1i*wc.Omega*t(j))));
                set(gca,'clim',[-0.5 0.5]);
                shading flat;
                axis equal;
                axis tight;
                title({'PntSrcWaveFieldUT - testHeaveCenter', '', 'Moving Elevation'});
                pause(0.01)
            end

            vel = wf.Velocity;

            inds = 1:10:201;
            X2 = X(inds, inds);
            Y2 = Y(inds, inds);

            U = vel{1}(inds, inds);
            V = vel{2}(inds, inds);

            subplot(2,1,2);
            quiver(X2,Y2,real(U), real(V));
            axis equal;
            axis tight;
            title('Velocity');
        end
        
        function testHeaveNonCenter(testCase)
            rho = 1000;

            x = -10:0.1:10;
            y = -10:0.1:10;

            [X, Y] = meshgrid(x, y);

            wc = PlaneWaves(1, 1, 0, 10);
            loc = [4 4];

            wf = PntSrcWaveField('Heave', loc, rho, wc, 1, X, Y);

            eta = squeeze(wf.Elevation);
            eta = eta{1};

            figure;
            subplot(2,1,1);
            pcolor(X,Y,real(eta));
            set(gca,'clim',[-0.5 0.5]);
            shading flat;
            axis equal;
            axis tight;
            title({'PntSrcWaveFieldUT - testHeaveNonCenter', '', 'Elevation'});

            vel = wf.Velocity;

            inds = 1:10:201;
            X2 = X(inds, inds);
            Y2 = Y(inds, inds);

            U = vel{1}(inds, inds);
            V = vel{2}(inds, inds);

            subplot(2,1,2);
            quiver(X2,Y2,real(U), real(V));
            axis equal;
            axis tight;
            title('Velocity');
        end
        
        function testMultiComp(testCase)
            rho = 1000;

            x = -10:0.1:10;
            y = -10:0.1:10;

            [X, Y] = meshgrid(x, y);

            loc = [0 0];

            T = [4, 3, 2, 1];
            nT = length(T);

            wcs = PlaneWaves(ones(nT, 1), T, zeros(nT, 1), 10);

            wf = PntSrcWaveField('Heave', loc, rho, wcs, 1, X, Y);

            eta = wf.Elevation;

            figure;
            for n = 1:nT
                subplot(2,2,n);
                pcolor(X,Y,real(eta{n}));
                shading flat;
                axis equal;
                axis tight
                if (n == 1)
                    title({'PntSrcWaveFieldUT - testMultiComp', '', ['Elevation, T = ' num2str(T(n))]});
                else
                    title(['Elevation, T = ' num2str(T(nT))]);
                end
            end
        end
        
        function testSurge(testCase)
            rho = 1000;

            x = -10:0.1:10;
            y = -10:0.1:10;

            [X, Y] = meshgrid(x, y);

            wc = PlaneWaves(1, 1, 0, 10);
            loc = [0 0];

            wf = PntSrcWaveField('Surge', loc, rho, wc, 1, X, Y);

            eta = wf.Elevation;
            eta = eta{1};

            t = 0:0.1:10;
            figure;

            subplot(2,1,1);
            for j = 1:length(t)
                cla(gca);
                pcolor(X,Y,real(eta*exp(1i*wc.Omega*t(j))));
                set(gca,'clim',[-0.5 0.5]);
                shading flat;
                axis equal;
                axis tight;    
                title({'PntSrcWaveFieldUT - testSurge', '', 'Moving Elevation'});
                pause(0.01)
            end

            vel = squeeze(wf.Velocity);

            inds = 1:10:201;
            X2 = X(inds, inds);
            Y2 = Y(inds, inds);

            U = vel{1}(inds, inds);
            V = vel{2}(inds, inds);

            subplot(2,1,2);
            quiver(X2,Y2,real(U), real(V));
            axis equal;
            axis tight;
            title('Velocity');
        end
        
        function testNoVel(testCase)
            rho = 1000;

            x = -10:0.5:10;
            y = -10:0.5:10;
            loc = [0 0];

            [X, Y] = meshgrid(x, y);

            wc = PlaneWaves(1, 2, 0, 10);

            % default with velocity
            wf = PntSrcWaveField('Surge', loc, rho, wc, 1, X, Y);
            
            testCase.verifyTrue(wf.HasVelocity);
            
            try
                vel = wf.Velocity;
                % Should not cause an error
            catch
                testCase.verifyFail();
            end
            
            pnts = [0 0 0; 1 1 0];
            % default 2 with velocity
            wf = PntSrcWaveField('Surge', loc, rho, wc, 0, pnts);
            
            testCase.verifyTrue(wf.HasVelocity);
            
            try
                vel = wf.Velocity;
                % Should not cause an error
            catch
                testCase.verifyFail();
            end
            
            % no velocity case
           wf = PntSrcWaveField('Surge', loc, rho, wc, 1, X, Y, 'NoVel');
            
            testCase.verifyTrue(~wf.HasVelocity);
            
            try
                vel = wf.Velocity;
                % Should cause an error
                testCase.verifyFail();
            catch
                % Pass - does throw error
            end
            
            wf = PntSrcWaveField('Surge', loc, rho, wc, 0, pnts, 'NoVel');
            
            testCase.verifyTrue(~wf.HasVelocity);
            
            try
                vel = wf.Velocity;
                testCase.verifyFail();
                % Should cause an error
            catch
                % Pass - does throw error
            end
        end
    end
end