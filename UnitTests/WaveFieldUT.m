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
classdef WaveFieldUT < matlab.unittest.TestCase
    % Unit test for WaveCompClass

    methods (Test)
       
        function test1(testCase)
            rho = 1000;

            x = -10:0.1:10;
            y = -10:0.1:10;

            [X Y] = meshgrid(x, y);

            T = [5, 4, 3, 2, 1];
            beta = pi/4;
            nT = length(T);

            wcs = PlaneWaves(ones(nT, 1), T, beta*ones(nT, 1), 10);

            wfin = PlaneWaveField(rho, wcs, 1, X, Y);
            wfb = PlaneWaveField(rho, wcs, 1, X, Y);
            
            testCase.verifyEqual(wfin, wfb);

            wc = PlaneWaves(1, 1, 0, 10);

            wfout = PlaneWaveField(rho, wc, 1, X, Y);
            
            testCase.verifyNotEqual(wfout, wfb);
        end
        
        function testSum(testCase)
            rho = 1000;

            x = -10:0.1:10;
            y = -10:0.1:10;

            [X, Y] = meshgrid(x, y);

            wc1 = PlaneWaves(1, 2, 0, 10);
            wc2 = PlaneWaves(0.5, 2, 0, 10);

            wfin = PlaneWaveField(rho, wc1, 1, X, Y);
            wfb = PlaneWaveField(rho, wc2, 1, X, Y);

            wfout = wfin + wfb;

            etaIn = wfin.Elevation;
            etaB = wfb.Elevation;
            etaOut = wfout.Elevation;

            figure;
            subplot(2,2,1);
            pcolor(X,Y,real(etaIn{1}));
            shading flat;
            axis equal;
            axis tight;
            set(gca, 'clim', [-1.5 1.5]);
            colorbar;

            subplot(2,2,2);
            pcolor(X,Y,real(etaB{1}));
            shading flat;
            axis equal;
            axis tight;
            set(gca, 'clim', [-1.5 1.5]);
            colorbar;

            subplot(2,2,3);
            pcolor(X,Y,real(etaOut{1}));
            shading flat;
            axis equal;
            axis tight;
            set(gca, 'clim', [-1.5 1.5]);
            colorbar;

            wc3 = PlaneWaves(1, 2, 0, 10);
            wc3.Epsilon = pi;

            wfd = PlaneWaveField(rho, wc3, 1, X, Y);

            wfe = wfin + wfd;

            etaD = wfd.Elevation;
            etaE = wfe.Elevation;

            figure;
            subplot(2,2,1);
            pcolor(X,Y,real(etaIn{1}));
            shading flat;
            axis equal;
            axis tight;
            set(gca, 'clim', [-1 1]);
            colorbar;

            subplot(2,2,2);
            pcolor(X,Y,real(etaD{1}));
            shading flat;
            axis equal;
            axis tight;
            set(gca, 'clim', [-1 1]);
            colorbar;

            subplot(2,2,3);
            pcolor(X,Y,real(etaE{1}));
            shading flat;
            axis equal;
            axis tight;
            set(gca, 'clim', [-1 1]);
            colorbar;
        end
        
        function testSpectrum(testCase)
            rho = 1000;
    
            x = -10:0.1:10;
            y = -10:0.1:10;

            [X, Y] = meshgrid(x, y);

            Hs = 2;
            fm = 0.2;
            df = 0.01;
            f = 0.1:df:1;
            S = bretschneider(Hs, fm, f).';

            a = sqrt(2*df*S);
            ep = 2*pi*rand(size(a));
            h = 10;
            
            wcs = PlaneWaves(a, 1./f, zeros(size(a)), h);

            wf = PlaneWaveField(rho, wcs, 1, X, Y);

            Spec = wf.Spectra('Points', [0 0]);

            figure;
            plot(Spec.Frequencies, Spec.Spectrum);
            hold on;
            plot(f, S, 'k--');
            
            testCase.verifyEqual(Spec.Spectrum, S, 'AbsTol', 1e-5);
        end
        
        function testSubtract(testCase)
            rho = 1000;

            x = -10:0.1:10;
            y = -10:0.1:10;

            [X, Y] = meshgrid(x, y);

            wc1 = PlaneWaves(1, 2, 0, 10);
            wc2 = PlaneWaves(1, 2, 0, 10);

            wfin = PlaneWaveField(rho, wc1, 1, X, Y);
            wfb = PlaneWaveField(rho, wc2, 1, X, Y);

            wfout = wfin - wfb;

            etaIn = wfin.Elevation;
            etaB = wfb.Elevation;
            etaOut = wfout.Elevation;

            figure;
            subplot(2,2,1);
            pcolor(X,Y,real(etaIn{1}));
            shading flat;
            axis equal;
            axis tight;
            set(gca, 'clim', [-1.5 1.5]);
            colorbar;

            subplot(2,2,2);
            pcolor(X,Y,real(etaB{1}));
            shading flat;
            axis equal;
            axis tight;
            set(gca, 'clim', [-1.5 1.5]);
            colorbar;

            subplot(2,2,3);
            pcolor(X,Y,real(etaOut{1}));
            shading flat;
            axis equal;
            axis tight;
            set(gca, 'clim', [-1.5 1.5]);
            colorbar;
            
            testCase.verifyEqual(real(etaOut{1}), zeros(size(X)), 'AbsTol', 1e-5);
        end
        
        function testMultiply(testCase)
            rho = 1000;

            x = -10:0.1:10;
            y = -10:0.1:10;

            [X, Y] = meshgrid(x, y);

            wc1 = PlaneWaves(1, 2, 0, 10);
  
            wfin = PlaneWaveField(rho, wc1, 1, X, Y);
            A = complex(1,1);

            wfout = A.*wfin;

            etaIn = wfin.Elevation;
            etaOut = wfout.Elevation;

            figure;
            subplot(1,2,1);
            pcolor(X,Y,real(etaIn{1}));
            shading flat;
            axis equal;
            axis tight;
            set(gca, 'clim', [-1.5 1.5]);
            colorbar;

            subplot(1,2,2);
            pcolor(X,Y,real(etaOut{1}));
            shading flat;
            axis equal;
            axis tight;
            set(gca, 'clim', [-1.5 1.5]);
            colorbar;
        end
        
        function testFlux(testCase)
            % compute flux across a unit width wave front
    
            % wave field
            rho = 1000;

            pnts = [0 0 0];
            wc1 = PlaneWaves(1, 2, 0, 10);

            wf = PlaneWaveField(rho, wc1, 0, pnts);

            % surface
            nrms = [1 0 0];
            ars = 1;
            surf = ControlSurface(pnts, nrms, ars);

            F = wf.EnergyFlux(surf);

            Fex = wc1.EnergyFlux(rho);
            
            testCase.verifyEqual(F{1,1}, Fex, 'AbsTol', 1e-4);

            % compute fluxes at line of 11 points at x = 0, y = -5:5
            pnts = zeros(11,3);
            pnts(:,2) = (-5:5)';

            wc2 = PlaneWaves([1 1], [2 4], [0 0], 10);

            wf = PlaneWaveField(rho, wc2, 0, pnts);

            nrms = zeros(11,3);
            nrms(:,1) = ones(11,1);
            ars = ones(11,1);

            surf = ControlSurface(pnts, nrms, ars);

            F = wf.EnergyFlux(surf, 'PerPoint');

            Fex = wc2.EnergyFlux(rho);
            
            testCase.verifyEqual(F{1}(1), Fex(1), 'AbsTol', 1e-4);
            testCase.verifyEqual(F{2}(1), Fex(2), 'AbsTol', 1e-4);

            % compute fluxes through a circle

            surf = CirContSurf([0, 0], 5, 50);

            wc = PlaneWaves(1, 2, 0, 10);

            wf = PlaneWaveField(rho, wc, 0, surf.Points);

            F = wf.EnergyFlux(surf);
            
            testCase.verifyEqual(F{1}, 0, 'AbsTol', 1e-6);

            % Flux though a circle with interpolation

            surf = CirContSurf([0, 0], 5, 50);

            wc = PlaneWaves(1, 2, 0, 10);

            x = -10:0.1:10;
            y = -10:0.1:10;

            [X Y] = meshgrid(x, y);

            wf = PlaneWaveField(rho, wc, 1, X, Y);

            F = wf.EnergyFlux(surf);
            
            testCase.verifyEqual(F{1}, 0, 'AbsTol', 1e-6);
        end
        
        function testRemove(testCase)
            rho = 1000;

            x = -10:0.1:10;
            y = -10:0.1:10;

            [X, Y] = meshgrid(x, y);

            wc = PlaneWaves(1, 1, 0, 10);
            loc = [0 0];

            wf1 = PntSrcWaveField('Heave', loc, rho, wc, 1, X, Y);
            wf2 = PntSrcWaveField('Heave', loc, rho, wc, 1, X, Y);
            
            eta1 = wf1.Elevation;
            
            thet = linspace(0,2*pi,21)';
            R = 2;
            cir = R*[cos(thet), sin(thet)];
            
            wf1.RemoveGeometries(cir, 'Out');
            eta2 = wf1.Elevation;
            
            wf2.RemoveGeometries(cir, 'In');
            eta3 = wf2.Elevation;
                        
            figure;
            subplot(3,1,1);
            pcolor(X,Y,real(eta1{1}))
            fet;
            set(gca,'clim', [-0.2 0.2]);
            title({'WaveFieldUT - testRemove', 'Full WaveField'});
            
            subplot(3,1,2);
            pcolor(X,Y,real(eta2{1}))
            fet;
            set(gca,'clim', [-0.2 0.2]);
            title('WaveField with removed circle');
            
            subplot(3,1,3);
            pcolor(X,Y,real(eta3{1}))
            fet;
            set(gca,'clim', [-0.2 0.2]);
            title('WaveField inside circle');
        end
    end
end