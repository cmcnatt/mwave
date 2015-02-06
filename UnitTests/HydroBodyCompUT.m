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
classdef HydroBodyCompUT < matlab.unittest.TestCase
    % Unit test for HydroBodyComp class

    methods (Test)
        function testDiffMatSize(testCase)
            % 1 T, 1 Beta, 1 dof
            load([mwavePath '\UnitTests\files\hbc_ut_1_1_1']);
            
            a = 1;
            t = hbc_ut_1_1_1_hydroF.T;
            beta = hbc_ut_1_1_1_hydroF.Beta;
            h = hbc_ut_1_1_1_hydroF.H;
            iwaves = PlaneWaves(a, t, beta, h);
            
            hb = HydroBodyComp(hbc_ut_1_1_1_hydroF, hbc_ut_1_1_1_floatB);
            
            testCase.verifyEqual(size(hb.T), [1 1]);
            testCase.verifyEqual(size(hb.IncWaves), [1 1]);
            testCase.verifyEqual(hb.DoF, 1);
            testCase.verifyEqual(size(hb.A), [1 1]);
            testCase.verifyEqual(size(hb.B), [1 1]);
            testCase.verifyEqual(size(hb.Fex), [1 1]);
            
            hb.SetDpto(1e3);
            X = hb.Motions;
            V = hb.Velocities;
            P = hb.Power;
            
            Xex =  -0.2802 - 0.7301i;
            omega = 2*pi./hb.T;
            Vex = 1i*omega.*Xex;
            Pex = 1.9314e+03;
            
            testCase.verifyEqual(X, Xex, 'AbsTol', 1e-4);
            testCase.verifyEqual(V, Vex, 'AbsTol', 1e-4);
            testCase.verifyEqual(P, Pex, 'AbsTol', 1);
            
            % 1 T, 1 Beta, Multiple dof
            load([mwavePath 'UnitTests\files\hbc_ut_1_1_dof']);
            hb = HydroBodyComp (hbc_ut_1_1_dof_hydroF, hbc_ut_1_1_dof_floatB);
            
            testCase.verifyEqual(size(hb.T), [1 1]);
            testCase.verifyEqual(size(hb.IncWaves), [1 1]);
            testCase.verifyEqual(hb.DoF, 6);
            testCase.verifyEqual(size(hb.A), [1 6 6]);
            testCase.verifyEqual(size(hb.B), [1 6 6]);
            testCase.verifyEqual(size(hb.Fex), [1 1 6]);
            
            D = zeros(6,6);
            D(3,3) = 1e3;
            hb.SetDpto(D);
            X = squeeze(hb.Motions);
            V = squeeze(hb.Velocities);
            P = squeeze(hb.Power);
            
            Xex = [0.0012 - 0.6174i; 0; -0.2802 - 0.7301i; 0; 0.0008 - 0.4471i; 0];
            omega = 2*pi./hb.T;
            Vex = 1i*omega.*Xex;
            Pex = 1.9314e+03;
            
            testCase.verifyEqual(X, Xex, 'AbsTol', 1e-4);
            testCase.verifyEqual(V, Vex, 'AbsTol', 1e-3);
            testCase.verifyEqual(P(3), Pex, 'AbsTol', 1);
            
            % 1 T, Multiple Beta, 1 dof
            load([mwavePath 'UnitTests\files\hbc_ut_1_b_1']);
            
            a = 1;
            t = hbc_ut_1_b_1_hydroF.T;
            h = hbc_ut_1_b_1_hydroF.H;
            
            beta = hbc_ut_1_b_1_hydroF.Beta;
            iwaves(length(beta), 1) = PlaneWaves;
            for n = 1:length(beta)
                iwaves(n) = PlaneWaves(a, t, beta(n), h);
            end
            
            hb = HydroBodyComp (hbc_ut_1_b_1_hydroF, hbc_ut_1_b_1_floatB);
            
            nBeta = 5;
            testCase.verifyEqual(size(hb.T), [1 1]);
            testCase.verifyEqual(size(hb.IncWaves), [1 nBeta]);
            testCase.verifyEqual(hb.DoF, 1);
            testCase.verifyEqual(size(hb.A), [1 1]);
            testCase.verifyEqual(size(hb.B), [1 1]);
            testCase.verifyEqual(size(hb.Fex), [1 nBeta]);
            
            D = 1e3;
            hb.SetDpto(D);
            X = squeeze(hb.Motions);
            V = squeeze(hb.Velocities);
            P = squeeze(hb.Power);
            
            Xex = (-0.2802 - 0.7301i)*ones(1, nBeta);
            omega = 2*pi./hb.T;
            Vex = 1i*omega.*Xex;
            Pex = 1.9314e+03*ones(1, nBeta);

            testCase.verifyEqual(X, Xex, 'AbsTol', 1e-4);
            testCase.verifyEqual(V, Vex, 'AbsTol', 1e-3);
            testCase.verifyEqual(P, Pex, 'AbsTol', 1);
            
            % 1 T Multiple Beta, Multi dof
            load([mwavePath 'UnitTests\files\hbc_ut_1_b_dof']);
            hb = HydroBodyComp (hbc_ut_1_b_dof_hydroF, hbc_ut_1_b_dof_floatB);
            
            nBeta = 5;
            testCase.verifyEqual(size(hb.T), [1 1]);
            testCase.verifyEqual(size(hb.IncWaves), [1 nBeta]);
            testCase.verifyEqual(hb.DoF, 6);
            testCase.verifyEqual(size(hb.A), [1 6 6]);
            testCase.verifyEqual(size(hb.B), [1 6 6]);
            testCase.verifyEqual(size(hb.Fex), [1 nBeta 6]);
            
            D = zeros(6,6);
            D(3,3) = 1e3;
            hb.SetDpto(D);
            X = squeeze(hb.Motions);
            V = squeeze(hb.Velocities);
            P = squeeze(hb.Power);
            
            Xex = [0.0012 - 0.6174i, 0, -0.2802 - 0.7301i, 0, 0.0008 - 0.4471i, 0;...
                0.0010 - 0.5347i, 0.0006 - 0.3087i, -0.2802 - 0.7301i, -0.0004 + 0.2235i, 0.0007 - 0.3872i, 0;...
                0.0008 - 0.4366i, 0.0008 - 0.4366i, -0.2802 - 0.7301i, -0.0005 + 0.3161i, 0.0005 - 0.3161i, 0;...
                0.0006 - 0.3087i, 0.0010 - 0.5347i, -0.2802 - 0.7301i, -0.0007 + 0.3872i, 0.0004 - 0.2235i, 0;...
                0, 0.0012 - 0.6174i, -0.2802 - 0.7301i, -0.0008 + 0.4471i, -0.0000 + 0.0000i, 0];

            omega = 2*pi./hb.T;
            Vex = 1i*omega.*Xex;
            Pex = 1.9314e+03;
            
            testCase.verifyEqual(X, Xex, 'AbsTol', 1e-4);
            testCase.verifyEqual(V, Vex, 'AbsTol', 1e-3);
            for n = 1:nBeta
                testCase.verifyEqual(P(n,3), Pex, 'AbsTol', 1);
            end
            
            % Multi T, 1 Beta, 1 dof
            load([mwavePath 'UnitTests\files\hbc_ut_T_1_1']);
            
            nT = length(hbc_ut_T_1_1_hydroF.T);
            a = ones(nT,1);
            t = hbc_ut_T_1_1_hydroF.T;
            beta = hbc_ut_T_1_1_hydroF.Beta*ones(nT,1);
            h = hbc_ut_T_1_1_hydroF.H;
            iwaves = PlaneWaves(a, t, beta, h);
            
            hb = HydroBodyComp (hbc_ut_T_1_1_hydroF, hbc_ut_T_1_1_floatB);
            
            nT = 9;
            testCase.verifyEqual(size(hb.T), [nT 1]);
            testCase.verifyEqual(size(hb.IncWaves), [1 1]);
            testCase.verifyEqual(hb.DoF, 1);
            testCase.verifyEqual(size(hb.A), [nT 1]);
            testCase.verifyEqual(size(hb.B), [nT 1]);
            testCase.verifyEqual(size(hb.Fex), [nT 1]);
            
            hb.SetDpto(1e3);
            X = hb.Motions;
            V = hb.Velocities;
            P = hb.Power;
            
            Xex =  [-0.0000 - 0.0001i; -0.0118 - 0.0093i; -0.1301 - 0.0971i;...
                -0.2802 - 0.7301i; 0.7327 - 1.0155i; 0.9929 - 0.5742i; ...
                1.0083 - 0.3811i; 1.0052 - 0.2874i; 1.0025 - 0.2330i];
            omega = 2*pi./hb.T;
            Vex = 1i*omega.*Xex;
            Pex = 1e3*[0; 0.002; 0.1301; 1.9314; 3.4393; 2.1199; 1.4334; 1.0655; 0.8364];
            
            testCase.verifyEqual(X, Xex, 'AbsTol', 1e-4);
            testCase.verifyEqual(V, Vex, 'AbsTol', 1e-3);
            testCase.verifyEqual(P, Pex, 'AbsTol', 1);
            
            % Multi T, 1 Beta, Multi dof
            load([mwavePath 'UnitTests\files\hbc_ut_T_1_dof']);
            hb = HydroBodyComp (hbc_ut_T_1_dof_hydroF, hbc_ut_T_1_dof_floatB);
            
            nT = 9;
            testCase.verifyEqual(size(hb.T), [nT 1]);
            testCase.verifyEqual(size(hb.IncWaves), [1 1]);
            testCase.verifyEqual(hb.DoF, 6);
            testCase.verifyEqual(size(hb.A), [nT 6 6]);
            testCase.verifyEqual(size(hb.B), [nT 6 6]);
            testCase.verifyEqual(size(hb.Fex), [nT 1 6]);
            
            D = zeros(6,6);
            D(3,3) = 1e3;
            hb.SetDpto(D);
            X = squeeze(hb.Motions);
            V = squeeze(hb.Velocities);
            P = squeeze(hb.Power);
            
            Xex =  [0.0263 - 0.0636i, 0,  -0.0000 - 0.0001i, 0, 0.0630 - 0.1525i, 0;...
                0.0071 - 0.2825i, 0, -0.0118 - 0.0093i, 0, 0.0103 - 0.4122i, 0;...
                0.0026 - 0.4795i, 0, -0.1301 - 0.0971i, 0, 0.0024 - 0.4615i, 0;...
                0.0012 - 0.6174i, 0, -0.2802 - 0.7301i, 0, 0.0008 - 0.4471i, 0;...
                0.0005 - 0.7106i, 0, 0.7327 - 1.0155i, 0, 0.0003 - 0.4338i, 0;...
                0.0002 - 0.7758i, 0, 0.9929 - 0.5742i, 0, 0.0001 - 0.4409i, 0;...
                0.0001 - 0.8274i, 0, 1.0083 - 0.3811i, 0, 0.0001 - 0.4830i, 0;...
                0.0001 - 0.8750i, 0, 1.0052 - 0.2874i, 0, 0.0000 - 0.5924i, 0;...
                0.0001 - 0.9202i, 0, 1.0025 - 0.2330i, 0, 0.0001 - 0.8841i, 0];
            omega = 2*pi./hb.T;
            Vex = (1i*omega*ones(1, 6)).*Xex;
            Pex = 1e3*[0; 0.002; 0.1301; 1.9314; 3.4393; 2.1199; 1.4334; 1.0655; 0.8364];
            
            testCase.verifyEqual(X, Xex, 'AbsTol', 1e-4);
            testCase.verifyEqual(V, Vex, 'AbsTol', 1e-3);
            testCase.verifyEqual(P(:,3), Pex, 'AbsTol', 1);
            
            % Multi T, Multi Beta, 1 dof
            load([mwavePath 'UnitTests\files\hbc_ut_T_b_1']);
            nB = length(hbc_ut_T_b_1_hydroF.Beta);
            nT = length(hbc_ut_T_b_1_hydroF.T);
            
            a = ones(nT, 1);
            t = hbc_ut_T_b_1_hydroF.T;
            h = hbc_ut_T_b_1_hydroF.H;
            
            beta = hbc_ut_T_b_1_hydroF.Beta;
            iwaves(nB, 1) = PlaneWaves;
            for n = 1:nB
                iwaves(n) = PlaneWaves(a, t, beta(n)*ones(nT, 1), h);
            end
            
            hb = HydroBodyComp (hbc_ut_T_b_1_hydroF, hbc_ut_T_b_1_floatB);
            
            nT = 9;
            nBeta = 5;
            testCase.verifyEqual(size(hb.T), [nT 1]);
            testCase.verifyEqual(size(hb.IncWaves), [1 nBeta]);
            testCase.verifyEqual(hb.DoF, 1);
            testCase.verifyEqual(size(hb.A), [nT 1]);
            testCase.verifyEqual(size(hb.B), [nT 1]);
            testCase.verifyEqual(size(hb.Fex), [nT nBeta]);

            D = 1e3;
            hb.SetDpto(D);
            X = squeeze(hb.Motions);
            V = squeeze(hb.Velocities);
            P = squeeze(hb.Power);
            
            Xex = [-0.0000 - 0.0001i; -0.0118 - 0.0093i; -0.1301 - 0.0971i;...
                -0.2802 - 0.7301i; 0.7327 - 1.0155i; 0.9929 - 0.5742i; ...
                1.0083 - 0.3811i; 1.0052 - 0.2874i; 1.0025 - 0.2330i]*ones(1, nBeta);
            omega = 2*pi./hb.T;
            Vex = (1i*omega*ones(1, nBeta)).*Xex;
            Pex = 1e3*[0; 0.002; 0.1301; 1.9314; 3.4393; 2.1199; 1.4334; 1.0655; 0.8364]*ones(1, nBeta);
            
            testCase.verifyEqual(X, Xex, 'AbsTol', 1e-4);
            testCase.verifyEqual(V, Vex, 'AbsTol', 1e-3);
            testCase.verifyEqual(P, Pex, 'AbsTol', 1);
            
            % Multi T, Mult Beta, Mult dof
            load([mwavePath 'UnitTests\files\hbc_ut_T_b_dof']);
            hb = HydroBodyComp(hbc_ut_T_b_dof_hydroF, hbc_ut_T_b_dof_floatB);
            
            nT = 9;
            nBeta = 5;
            testCase.verifyEqual(size(hb.T), [nT 1]);
            testCase.verifyEqual(size(hb.IncWaves), [1 nBeta]);
            testCase.verifyEqual(hb.DoF, 6);
            testCase.verifyEqual(size(hb.A), [nT 6 6]);
            testCase.verifyEqual(size(hb.B), [nT 6 6]);
            testCase.verifyEqual(size(hb.Fex), [nT nBeta 6]);

            D = zeros(6, 6);
            D(3,3) = 1e3;
            hb.SetDpto(D);
            X = squeeze(hb.Motions);
            V = squeeze(hb.Velocities);
            P = squeeze(hb.Power);
            
            Pex = 1e3*[0; 0.002; 0.1301; 1.9314; 3.4393; 2.1199; 1.4334; 1.0655; 0.8364]*ones(1, nBeta);
            
            testCase.verifyEqual(size(X), [nT nBeta 6]);
            testCase.verifyEqual(size(V), [nT nBeta 6]);
            testCase.verifyEqual(squeeze(P(:,:,3)), Pex, 'AbsTol', 1);
        end
        
        function test1CylComp(testCase)
            % Test the array comptutation on a 2 cylinder array
            
            load([mwavePath 'UnitTests\files\cyl6dof_hbUT']);
            load([mwavePath 'UnitTests\files\cyl1_verUT']);

            wComp = wComp1cyl;
            waveW = wWave1cyl;
            floatB = wComp1cyl.Bodies;

            waveW.BodyMotions = wComp.Motions;

            hb = HydroBody(cylHydBod);
            hb.XYpos = floatB.XYpos;
            
            iwaves = wComp.IncWaves;
            bComp = HydroBodyComp(hb, iwaves);
            
            % Round the values to get rid of erroneous small values
            pwr10 = -1;
            
            A = round2val(squeeze(bComp.A), pwr10);
            B = round2val(squeeze(bComp.B), pwr10);
            Fex = round2val(squeeze(bComp.Fex), pwr10);

            Aw = round2val(squeeze(wComp.A), pwr10);
            Bw = round2val(squeeze(wComp.B), pwr10);
            Fexw = round2val(squeeze(wComp.Fex), pwr10);

%             errA = matErr(A, Aw);
%             errB = matErr(B, Bw);
%             errFex = matErr(Fex, Fexw);

            testCase.verifyEqual(size(A), [6 6]);
            testCase.verifyEqual(size(B), [6 6]);
            testCase.verifyEqual(size(Fex), [5 6]);
            
            for m = 1:6
                for n = 1:6
                    if (m ~= 6)
                        testCase.verifyEqual(Fex(m,n), Fexw(m,n), 'RelTol', 10^-2);         % less than 1% error
                    end
                    testCase.verifyEqual(A(m,n), Aw(m,n), 'RelTol', 3*10^-2);   % less than 3% error
                    testCase.verifyEqual(B(m,n), Bw(m,n), 'RelTol', 10^-2);     % less than 1% error
                end
            end
            
            [X, Y] = waveW.FieldPoints;
            waveB = bComp.WaveField(1, X, Y);
            
            etaW = waveW.Elevation('Total');
            etaB = waveB.Elevation('Total');
            
            figure;
            subplot(2,1,1)
            pcolor(X,Y,abs(etaW{1}));
            fet;
            set(gca, 'clim', [0.8 1.2], 'xtick', [], 'ytick', []);
            title({'HydroBodyCompUT', 'Wamit Total Wave Field'});
            
            subplot(2,1,2)
            pcolor(X,Y,abs(etaB{1}));
            fet;
            set(gca, 'clim', [0.8 1.2], 'xtick', [], 'ytick', []);
            title({'HydroBodyComp Total Wave Field'});
        end
    end
end