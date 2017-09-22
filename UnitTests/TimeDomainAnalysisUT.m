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
classdef TimeDomainAnalysisUT < matlab.unittest.TestCase

    methods (Test)
       
        function testSetGet(testCase)

            endTime = 100;
            dt = 0.01;
            Nsamp = endTime/dt;
            freqs = 1./(0.8:0.1:1.6);
            Nf = length(freqs);
            dof = 6;
            
            mots = zeros(Nf, dof, Nsamp);
            
            for m = 1:Nf
                for n = 1:dof
                    [mots(m, n, :), timeMot] = TimeDomainAnalysisUT.GenerateRandSig(freqs(m), dt, Nsamp);
                end
            end
            
            dt = 0.005;
            Nsamp = endTime/dt;
            ptoKin = zeros(Nf, 1, Nsamp);
            ptoDyn = zeros(Nf, 1, Nsamp);
            for m = 1:Nf
                [ptoKin(m, 1, :), timePto] = TimeDomainAnalysisUT.GenerateRandSig(freqs(m), dt, Nsamp);
                [ptoDyn(m, 1, :)] = TimeDomainAnalysisUT.GenerateRandSig(freqs(m), dt, Nsamp);
            end
            
            wgPos = [-2, 0];
            dt = 0.005;
            Nsamp = endTime/dt;
            waves = zeros(Nf, 1, Nsamp);
            for m = 1:Nf
                [waves(m, 1, :), timeWave] = TimeDomainAnalysisUT.GenerateRandSig(freqs(m), dt, Nsamp);
            end
                       
            tda = TimeDomainAnalysis;
            
            tda.SetMotions(1:6, timeMot, mots);
            tda.SetPtoKinematic(1, timePto, ptoKin);
            tda.SetPtoDynamic(1, timePto, ptoDyn);
            pow = ptoDyn.*ptoKin;
            tda.SetWaves(wgPos, timeWave, waves);
            
            motsOut = tda.GetMotions(1:Nf, 1:6);
            ptoKinOut = tda.GetPtoKinematic(1:Nf, 1);
            ptoDynOut = tda.GetPtoDynamic(1:Nf, 1);
            powOut = tda.Power(1:Nf, 1);
            wavesOut = tda.GetWaves(1:Nf, 1);
            
            for m = 1:Nf
                testCase.assertEqual(squeeze(ptoKin(m, 1, :)), ptoKinOut{m, 1});
                testCase.assertEqual(squeeze(ptoDyn(m, 1, :)), ptoDynOut{m, 1});
                testCase.assertEqual(squeeze(pow(m, 1, :)), powOut{m, 1});
                
                ptoKinOutm = tda.GetPtoKinematic(m, 1);
                ptoDynOutm = tda.GetPtoDynamic(m, 1);
                testCase.assertEqual(squeeze(ptoKin(m, 1, :)), ptoKinOutm{1});
                testCase.assertEqual(squeeze(ptoDyn(m, 1, :)), ptoDynOutm{1});
                
                testCase.assertEqual(squeeze(waves(m, 1, :)), wavesOut{m, 1});
                wavesOutm = tda. GetWaves(m, 1);
                testCase.assertEqual(squeeze(waves(m, 1, :)), wavesOutm{1});
                
                for n = 1:dof
                    testCase.assertEqual(squeeze(mots(m, n, :)), motsOut{m, n});
                    motsOutmn = tda.GetMotions(m, n);
                    testCase.assertEqual(squeeze(mots(m, n, :)), motsOutmn{1});
                end
            end
        end
        
        function testWaveAmps(testCase)
            wgPos = [-2, 0];
            beta = [0 pi];
            Nbeta = length(beta);
            
            Nf = 11;
            waveA = zeros(Nf, 1, Nbeta);
            
            for m = 1:Nf
                for n = 1:Nbeta
                    waveA(m, n) = abs(0.1*randn + 1)*exp(1i*2*pi*rand);
                end
            end
                       
            tda = TimeDomainAnalysis;
            
            tda.SetWaves(wgPos, [], waveA, beta);
            
            [waveAOut, ~, wgPosOut, betaOut] = tda.GetWaves(1:Nf, 1);
            
            testCase.assertEqual(wgPos, wgPosOut);
            testCase.assertEqual(beta, betaOut);
            
            for m = 1:Nf
                for n = 1:Nbeta
                    testCase.assertEqual(squeeze(waveA(m, 1, n)), waveAOut{m}(n));
                end
            end
        end
        
        function testSpectra(testCase)
            endTime = 100;
            dt = 0.01;
            Nsamp = endTime/dt;
            freqs = 1./(0.8:0.2:1.6);
            Nf = length(freqs);
            dof = 3;
            
            mots = zeros(Nf, dof, Nsamp);
            
            for m = 1:Nf
                for n = 1:dof
                    [mots(m, n, :), timeMot] = TimeDomainAnalysisUT.GenerateRandSig(freqs(m), dt, Nsamp);
                end
            end
            
            tda = TimeDomainAnalysis;
            
            tda.SetMotions(1:dof, timeMot, mots);
            
            [specMots, f] = tda.GetMotions(1:Nf, 1:3, 'spectra');
            
            for n = 1:dof
                figure;
                for m = 1:Nf
                    subplot(Nf,2,(m-1)*2+1)
                    plot(timeMot, squeeze(mots(m, n, :)));
                    subplot(Nf, 2, (m-1)*2+2)
                    plot(f{m, n}, abs(specMots{m, n}).^2);
                    set(gca, 'xlim', [0 6]);
                end
            end
            
        end
                
        function testTimeLims(testCase)
            endTime = 100;
            dt = 0.01;
            rampTime = 10;
            Nsamp = endTime/dt;
            freqs = 1./(0.8:0.2:1.6);
            Nf = length(freqs);
            dof = 3;
            
            mots = zeros(Nf, dof, Nsamp);
            
            for m = 1:Nf
                for n = 1:dof
                    [mots(m, n, :), timeMot] = TimeDomainAnalysisUT.GenerateRandSig(freqs(m), dt, Nsamp, 'rampTime', rampTime);
                end
            end
            
            tda = TimeDomainAnalysis;
            
            tda.SetMotions(1:dof, timeMot, mots);
            
            tda.SetMotionTimeLimits(1:Nf, 1:dof, 10, Inf);
            
            [motsOut, timeOut] = tda.GetMotions(1:Nf, 1:3);
            
            for n = 1:dof
                figure;
                for m = 1:Nf
                    subplot(Nf,2,(m-1)*2+1)
                    plot(timeMot, squeeze(mots(m, n, :)));
                    set(gca, 'xlim', [0 100]);
                    subplot(Nf, 2, (m-1)*2+2)
                    plot(timeOut{m, n}, motsOut{m, n});
                    set(gca, 'xlim', [0 100]);
                end
            end    
        end
        
        function testTimeLims2(testCase)
            endTime = 100;
            dt = 0.01;
            rampTime = 30;
            Nsamp = endTime/dt;
            freqs = 1./(0.8:0.2:1.6);
            Nf = length(freqs);
            dof = 3;
            
            mots = zeros(Nf, dof, Nsamp);
            
            for m = 1:Nf
                for n = 1:dof
                    [mots(m, n, :), timeMot] = TimeDomainAnalysisUT.GenerateRandSig(freqs(m), dt, Nsamp, 'rampTime', rampTime);
                end
            end
            
            tda = TimeDomainAnalysis;
            
            tda.SetMotions(1:dof, timeMot, mots);
            
            [specs1, f1] = tda.GetMotions(1:Nf, 1:3, 'spectra');
            
            tda.SetMotionTimeLimits(1:Nf, 1:dof, rampTime, Inf);
            [specs2, f2] = tda.GetMotions(1:Nf, 1:3, 'spectra');
            
            for n = 1:dof
                figure;
                for m = 1:Nf
                    subplot(Nf,2,(m-1)*2+1)
                    plot(f1{m,n}, abs(specs1{m,n}).^2);
                    set(gca, 'xlim', [0 5]);
                    subplot(Nf, 2, (m-1)*2+2)
                    plot(f2{m,n}, abs(specs2{m,n}).^2);
                    set(gca, 'xlim', [0 5]);
                end
            end
        end
        
        function testTimeLims3(testCase)
            endTime = 100;
            dt = 0.01;
            rampTime = 30;
            Nsamp = endTime/dt;
            freqs = 1./(0.8:0.2:1.6);
            Nf = length(freqs);
            dof = 1;
            
            ptoKin = zeros(Nf, dof, Nsamp);
            
            for m = 1:Nf
                for n = 1:dof
                    [ptoKin(m, n, :), timePto] = TimeDomainAnalysisUT.GenerateRandSig(freqs(m), dt, Nsamp, 'rampTime', rampTime);
                end
            end
            
            tda = TimeDomainAnalysis;
            
            tda.SetPtoKinematic(1:dof, timePto, ptoKin);
            tda.SetPtoDynamic(1:dof, timePto, ptoKin);
            
            [pow1, time1] = tda.Power(1:Nf, 1);
            
            tda.SetPtoTimeLimits(1:Nf, 1:dof, rampTime, Inf);
            [pow2, time2] = tda.Power(1:Nf, 1);
            
            for n = 1:dof
                figure;
                for m = 1:Nf
                    subplot(Nf,2,(m-1)*2+1)
                    plot(time1{m,n}, pow1{m,n});
                    set(gca, 'xlim', [0 100]);
                    subplot(Nf, 2, (m-1)*2+2)
                    plot(time2{m,n}, pow2{m,n});
                    set(gca, 'xlim', [0 100]);
                end
            end
        end
        
        function testRAO1(testCase)
            endTime = 100;
            dt = 0.01;
            rampTime = 30;
            Nsamp = endTime/dt;
            freqs = 1./(0.8:0.1:1.6);
            Nf = length(freqs);
            dof = 6;
            
            mots = zeros(Nf, dof, Nsamp);
            
            for m = 1:Nf
                for n = 1:dof
                    [mots(m, n, :), timeMot] = TimeDomainAnalysisUT.GenerateRandSig(freqs(m), dt, Nsamp, 'rampTime', rampTime);
                end
            end
            
            tda = TimeDomainRAO;
            tda.Frequencies = freqs;
            
            tda.SetMotions(1:dof, timeMot, mots);
            
            tda.SetMotionTimeLimits(1:Nf, 1:dof, rampTime, Inf);
            [raos, f] = tda.GetMotions(1:Nf, 1:dof, 'rao', 2);
            
            figure;
            for n = 1:dof
                subplot(dof,2,(n-1)*2+1)
                plot(f, abs(raos{n,1}), 'o-');
                
                subplot(dof, 2, (n-1)*2+2)
                plot(f, abs(raos{n,2}), 'o-');
            end
        end
        
        function testRAO2(testCase)
            endTime = 100;
            dt = 0.01;
            rampTime = 30;
            Nsamp = endTime/dt;
            freqs = 1./(0.8:0.1:1.6);
            Nf = length(freqs);
            dof = 1;
            
            ptoKin = zeros(Nf, dof, Nsamp);
            
            for m = 1:Nf
                for n = 1:dof
                    [ptoKin(m, n, :), timePto] = TimeDomainAnalysisUT.GenerateRandSig(freqs(m), dt, Nsamp, 'rampTime', rampTime);
                end
            end
            
            tda = TimeDomainRAO;
            tda.Frequencies = freqs;
            
            tda.SetPtoKinematic(1:dof, timePto, ptoKin);
            tda.SetPtoDynamic(1:dof, timePto, ptoKin);
            
            [pow1, f1] = tda.Power(1:Nf, 1, 'rao');
            
            tda.SetPtoTimeLimits(1:Nf, 1:dof, rampTime, Inf);
            [pow2, f2] = tda.Power(1:Nf, 1, 'rao');
            
            figure;
            for n = 1:dof
                subplot(dof,1,n)
                plot(f1, pow1{n}, 'o-');
                
                hold on
                plot(f2, pow2{n}, 'o-');
            end
        end
    end
    
    methods (Static)
        function [sig, time] = GenerateRandSig(freq, dt, Nsamp, varargin)
            
            [opts, args] = checkOptions({{'rampTime', 1}}, varargin);
            rampTime = [];
            if opts(1)
                rampTime = args{1};
            end
            
            time = dt:dt:(dt*Nsamp);
            
            Nf = length(freq);
            
            sig = zeros(size(time));
            for m = 1:Nf
                
                a(1) = abs(5*randn + 10);
                a(2) = abs(0.5*randn + 1);
                a(3) = abs(0.5*randn + 1);
                a(4) = abs(0.05*randn + 0.1);
                a(5) = abs(0.05*randn + 0.1);
                                
                for n = 1:5
                    ep = 2*pi*rand;
                    sig = sig + a(n)*cos(2*pi*n*freq(m)*time + ep);
                end
            end
            
            nA = abs(0.25*randn + 0.5);
            sig = sig + nA*rand(size(sig));
            if ~isempty(rampTime)
                ramp = ones(1, length(time));
                irampT = indexOf(time, rampTime);
                slope = 1/rampTime;
                ramp(1:irampT) = slope*time(1:irampT);
                sig = ramp.*sig;
            end
        end
    end
end