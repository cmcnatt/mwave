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
classdef WaveSpectrumUT < matlab.unittest.TestCase

    methods (Test)
       
        function testDirectional(testCase)
            Hs = 2;
            fm = 0.2;
            f = 0.1:0.01:1;
            Spec = Bretschneider(Hs, 1./fm, 1./f);
            S = Spec.Spectrum;

            s = 4;
            thetac = 0;
            theta = linspace(-pi, pi, 91);
            G = cosSpectSpread(s, thetac, theta);

            SG = S'*G; 
            Spec = WaveSpectrum(SG, f, theta);

            if (~Spec.IsDir)
                disp('FAIL: Spectrum should be directional but is not');
            end

            if (abs(Spec.SigWaveHeight - Hs) > 0.01)
                disp('FAIL: Spectrum Hs is not correct');
            end

            figure;
            pcolor(f, theta, Spec.Spectrum.');

            figure;
            plot(f, Spec.Spectrum('Nondir'));

            figure;
            pcolor(Spec.Frequencies('WithEnergy'), Spec.Directions('WithEnergy'), Spec.Spectrum('WithEnergy').');

            figure;
            pcolor(f, theta, Spec.Amplitudes.');
        end
        
        function testNondirectional(testCase)
            Hs = 2;
            fm = 0.2;
            f = 0.1:0.01:1;
            Spec = Bretschneider(Hs, 1./fm, 1./f);

            if (Spec.IsDir)
                disp('FAIL: Spectrum should not be directional but is.');
            end

            if (abs(Spec.SigWaveHeight - Hs) > 0.01)
                disp('FAIL: Spectrum Hs is not correct');
            end

            figure;
            plot(f, Spec.Spectrum);

            figure;
            plot(Spec.Frequencies('WithEnergy'), Spec.Spectrum('WithEnergy'));

            figure;
            plot(f, Spec.Amplitudes);
        end
    end
end