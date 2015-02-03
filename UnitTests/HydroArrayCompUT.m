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
classdef HydroArrayCompUT < matlab.unittest.TestCase
    % Unit test for HydroArray Class

    methods (Test)
                
        function test2CylComp(testCase)
            % Test the array comptutation on a 2 cylinder array
            
            load([mwavePath 'UnitTests\files\cyl6dof_hbUT']);
            load([mwavePath 'UnitTests\files\cyl2_verUT']);

            wComp = wComp2cyl;
            waveW = wWave2cyl;
            floatBs = wComp2cyl.Bodies;

            waveW.BodyMotions = wComp.Motions;

            for n = 1:length(floatBs)
                hbs(n) = HydroBody(cylHydBod);
                hbs(n).XYpos = floatBs(n).XYpos;
            end
            
            iwaves = wComp.IncWaves;
            aComp = HydroArrayComp(hbs, iwaves);
            
            % Round the values to get rid of erroneous small values
            pwr10 = 10;
            
            A = round2val(squeeze(aComp.A), pwr10);
            B = round2val(squeeze(aComp.B), pwr10);
            Fex = round2val(squeeze(aComp.Fex), pwr10);

            Aw = round2val(squeeze(wComp.A), pwr10);
            Bw = round2val(squeeze(wComp.B), pwr10);
            Fexw = round2val(squeeze(wComp.Fex), pwr10);

%             errA = matErr(A, Aw);
%             errB = matErr(B, Bw);
%             errFex = matErr(Fex, Fexw);

            testCase.verifyEqual(size(A), [12 12]);
            testCase.verifyEqual(size(B), [12 12]);
            testCase.verifyEqual(size(Fex), [5 12]);
            
            for m = 1:12
                for n = 1:5
                    testCase.verifyEqual(Fex(n,m), Fexw(n,m), 'RelTol', 10^-2);         % less than 1% error
                end
                for n = 1:12
                    testCase.verifyEqual(A(m,n), Aw(m,n), 'RelTol', 10^-2);   % less than 4% error (error only significant for relatively small values)
                    testCase.verifyEqual(B(m,n), Bw(m,n), 'RelTol', 10^-2);     % less than 1% error
                end
            end
            
            [X, Y] = waveW.FieldPoints;
            waveA = aComp.WaveField(1, X, Y);
            
            etaW = waveW.Elevation('Total');
            etaA = waveA.Elevation('Total');
            
            bet = {'0', '\pi/8', '\pi/4', '3\pi/8', '\pi/2'};
            figure;
            for n = 1:5
                subplot(5,2,(n-1)*2+1)
                pcolor(X,Y,abs(etaW{n}));
                fet;
                set(gca, 'clim', [0.8 1.2], 'xtick', [], 'ytick', []);
                if (n == 1)
                    title({'HydroArrayCompUT', 'Wamit Total Wave Field'});
                end
                ylabel(['\beta = ' bet{n}]);

                subplot(5,2,(n-1)*2+2)
                pcolor(X,Y,abs(etaA{n}));
                fet;
                set(gca, 'clim', [0.8 1.2], 'xtick', [], 'ytick', []);
                if (n == 1)
                    title({'HydroArrayComp Total Wave Field'});
                end
            end
        end
        
        function testIncCirWave(testCase)
            % Test the array comptutation on a 2 cylinder array
            
            load([mwavePath 'UnitTests\files\cyl6dof_hbUT']);

            t = cylHydBod.T;
            h = cylHydBod.H;
            pos = [0 -4; 0 4];
            for n = 1:2
                hbs(n) = HydroBody(cylHydBod);
                hbs(n).XYpos = pos(n,:);
            end
            
            s = 10;
            R = 2;
            randPhase = true;
            coefs = HydroArrayCompUT.cirCoefs(t, h, s, R, randPhase);
            iwaves = CirWaves('In', [0 0], coefs, t, h);
            aComp = HydroArrayComp(hbs, iwaves);
            
            x = -20:0.2:20;
            [X, Y] = meshgrid(x, x);
            
            waveA = aComp.WaveField(1, X, Y);
            etaA = waveA.Elevation('Total');
            
            figure;

            pcolor(X,Y,real(etaA{1}));
            fet;
            set(gca, 'clim', [-1.2 1.2], 'xtick', [], 'ytick', []);
            title({'HydroArrayCompUT', 'Short-crested Wave Field'});
        end
    end
    
    methods (Static)
        
        function [coefs] = cirCoefs(t, h, s, R, randPhase)
            
            coefs = cell(length(t),1);
            
            for n = 1:length(t)
                % wave component
                k = IWaves.SolveForK(2*pi./t(n), h);
                
                % directional coefficients
                Nthet = 4096;
                dthet = 2*pi/Nthet;
                theta = 0:dthet:(2*pi-dthet);

                %A = gamma(s + 1)/(gamma(s + 1/2)*2*sqrt(pi));
                D = cos(1/2*theta).^(2*s);
                A = trapz(theta,D);
                D = D./A;

%                 intD = trapz(theta, sqrt(D));
%                 R = 1/(k*(2*pi)^3*intD^2);
%                 R = 2;

                if (R < 0)
                    intD = trapz(theta, sqrt(D));
                    R = 1/(k*(2*pi)^3*intD^2);
                end
                F = (2*pi)^(3/2)*sqrt(k*R*D);
                if (randPhase)
                    F = F.*exp(1i*2*pi*rand(size(theta)));
                end

                fori = 1/Nthet*fft(F);
                M = 80;
                fm = ones(2*M+1,1);
                fm(M+1) = fori(1);
                fm(M+2:2*M+1) = fori(2:M+1);
                fm(M:-1:1) = fori(Nthet:-1:Nthet-M+1);
                fm = 2*pi*fm;
                m = (-M:M).';
                fm = (-1i).^m.*fm;
                
                coefs{n} = fm.';
            end
        end

    end
end

