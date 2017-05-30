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
classdef WaveFieldCollectionUT < matlab.unittest.TestCase

    methods (Test)
       
        function test1(testCase)
            % Create a directional wave field with a few periods and a
            % couple directions
            
            rho = 1000;

            x = -10:0.1:10;
            y = -10:0.1:10;

            [X, Y] = meshgrid(x, y);

            T = [5, 4, 3, 2, 1];
            beta = [ 0, pi/4];
            nT = length(T);
            nB = length(beta);

            for m = 1:nB
                wcs = PlaneWaves(ones(size(T)), T, beta(m)*ones(size(T)), 10);
                wfs(m) = PlaneWaveField(rho, wcs, 1, X, Y);
            end

            wf = WaveFieldCollection(wfs, 'Direction', beta);
            
            testCase.verifyEqual(wf.CollType, 'Direction');
            testCase.verifyEqual(wf.WFcount, nB);

            eta = wf.Elevation;

            figure;
            for m = 1:nB
                testCase.verifyEqual(wf.Indices(m), beta(m));
                for n = 1:nT
                    subplot(2,5,(m-1)*nT + n);
                    pcolor(X,Y,real(eta{n, m}));
                    if(m == 1)
                        if (n == 3)
                            title({'WaveFieldCollectionUT - test1', ['T = ' num2str(T(n))]});
                        else
                            title(['T = ' num2str(T(n))]);
                        end
                    end
                    if (n == 1)
                        ylabel(['Dir = ' num2str(beta(m))]);
                    end
                    set(gca, 'xtick', [], 'ytick', []);
                    shading flat;
                    axis equal;
                    axis tight
                end
            end
            
            vel = wf.Velocity;
            testCase.verifyEqual(size(vel), [nT, nB, 3]);

            hs = wf.SigWaveHeight('Merge');
            testCase.verifyEqual(size(hs), size(X));
            
            hs = wf.SigWaveHeight;
            testCase.verifyEqual(length(hs), nB);
        end
        
        function test2(testCase)
            rho = 1000;

            x = -10:0.1:10;
            y = -10:0.1:10;

            [X, Y] = meshgrid(x, y);

            Hs = 2;
            fm = 0.2;
            df = 0.1;
            f = 0.1:df:1;
            Spec = Bretschneider(Hs, 1./fm, 1./f);
            S = Spec.Spectrum;

            s = 10;
            thetac = pi/180*0;
            dtheta = 10;
            theta = pi/180*(-40:dtheta:40);
            G = cosSpectSpread(s, thetac, theta);

            SG = G'*S; 

            a = sqrt(2*df*dtheta*SG);
            ep = 2*pi*rand(size(a));
            h = 10;

            for n = 1:length(theta)
                wcs = PlaneWaves(a(n,:), 1./f, theta(n)*ones(size(f)), h, ep(n,:));
                wfs(n) = PlaneWaveField(rho, wcs, 1, X, Y);
            end

            wf = WaveFieldCollection(wfs, 'Direction', theta);


            Spec = wf.Spectra('Points', [0 0], 'Merge');

            figure;
            subplot(2,1,1)
            pcolor(Spec.Frequencies, Spec.Directions, Spec.Spectrum.');
            title({'WaveFieldCollectionUT - test2', 'WaveFieldCollection.Spectra(''Merge'')'});
            xlabel('Frequency (Hz)');
            ylabel('Directions (deg)');
            
            subplot(2,1,2)
            pcolor(f, theta, SG);
            title('Original Spectrum');
            xlabel('Frequency (Hz)');
            ylabel('Directions (deg)');
            
            SpecNoMerge = wf.Spectra('Points', [0 0]);
            
            Ntheta = length(theta);
            testCase.verifyEqual(length(SpecNoMerge), Ntheta);
            
            figure;
            for n = 1:Ntheta
                subplot(Ntheta, 1, n);
                plot(f, SpecNoMerge{n}.Spectrum);
                if (n == 1)
                    title({'WaveFieldCollectionUT - test2', 'WaveFieldCollection.Spectra'});
                end
                if (n == Ntheta)
                    xlabel('Frequency (Hz)');
                else
                    set(gca, 'xtick', []);
                end
            end
        end
        
        function test3(testCase)
            rho = 1000;

            x = -10:0.1:10;
            y = -10:0.1:10;

            [X, Y] = meshgrid(x, y);

            Hs = 2;
            fm = 0.2;
            df = 0.1;
            f = 0.1:df:1;
            Spec = Bretschneider(Hs, 1./fm, 1./f);
            S = Spec.Spectrum;

            s = 10;
            thetac = pi/180*0;
            dtheta = 10;
            theta = pi/180*(-40:dtheta:40);
            G = cosSpectSpread(s, thetac, theta);

            SG = G'*S; 

            a = sqrt(2*df*dtheta*SG);
            ep = 2*pi*rand(size(a));
            h = 10;

            for n = 1:length(theta)
                wcs = PlaneWaves(a(n,:), 1./f, theta(n)*ones(size(f)), h, ep(n,:));
                wfs(n) = PlaneWaveField(rho, wcs, 1, X, Y);
            end

            wfa = WaveFieldCollection(wfs, 'Direction', theta);
            wfb = WaveFieldCollection(wfs, 'Direction', theta);
            
            testCase.verifyEqual(wfa, wfb);
        end
        
        function test4(testCase)
            rho = 1000;
            x = -10:0.1:10;
            y = -10:0.1:10;
            h = 10;
            a = 1;
            T = [1 2 3 4];

            [X, Y] = meshgrid(x, y);
            
            wcs = PlaneWaves(ones(size(T)), T, zeros(size(T)), h);
            
            % surge 
            swave = PntSrcWaveField('Surge', [0 0], rho, wcs, 1, X, Y);
            % heave
            hwave = PntSrcWaveField('Heave', [0 0], rho, wcs, 1, X, Y);
            
            radwaves = WaveFieldCollection([swave hwave], 'Radiated', {'Surge', 'Heave'});
            
            beta = [0 pi/4 pi/2];
            
            a = 0.4;
            b = 0.6;
            motions = a + (b-a).*rand(length(T), length(beta), 2);
            rwfs(1, length(beta)) = WaveFieldCollection;
            
            for n = 1:length(beta)
                wcs.Beta = beta(n)*ones(wcs.Count,1);
                wfs(n) = PlaneWaveField(rho, wcs, 1, X, Y);
                rwfs(n) = radwaves;
            end
            
            incwaves = WaveFieldCollection(wfs, 'Direction', beta);
            radwavesM = WaveFieldCollection(rwfs, 'Direction', beta);
            radwavesM = motions*radwavesM;
            
            siz = radwavesM.GetTotalSize;
            sizExp = [length(T), length(beta), 2];
            
            testCase.verifyEqual(sizExp, siz);
            totwaves = incwaves + radwavesM;

            eta = totwaves.Elevation;
            
            figure;
            for m = 1:length(T)
                for n = 1:length(beta)
                    ind = (m-1)*length(beta)+n;
                    subplot(length(T), length(beta), ind)
                    pcolor(abs(eta{m,n,1}));
                    fet;
                    set(gca, 'clim', [0.8 1.2], 'xtick', [], 'ytick', []);
                    
                    if (ind == 2)
                        title({'WaveFieldCollectionUT - test4', 'Surge'})
                    end
                end
            end
            
            figure;
            for m = 1:length(T)
                for n = 1:length(beta)
                    ind = (m-1)*length(beta)+n;
                    subplot(length(T), length(beta), ind)
                    pcolor(abs(eta{m,n,2}));
                    fet;
                    set(gca, 'clim', [0.8 1.2], 'xtick', [], 'ytick', []);
                    
                    if (ind == 2)
                        title({'WaveFieldCollectionUT - test4', 'Heave'})
                    end
                end
            end
            
            vel = totwaves.Velocity;

            testCase.verifyEqual(size(vel), [length(T), length(beta), 2, 3]);
        end
    end
end