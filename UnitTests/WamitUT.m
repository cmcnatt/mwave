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
classdef WamitUT < matlab.unittest.TestCase

    methods (Test)
       
        function testRAO(testCase)
            name = 'ut1';
            folder = [WamitUT.getMainPath name];
            T = 0.1:0.1:6;

            cyl = WamitUT.createCylFB();
            cyl.Handle = 'RAO';

            ut1_run = WamitRunCondition(folder, name);
            ut1_run.Rho = WamitUT.getRho;

            ut1_run.FloatingBodies = cyl;
            ut1_run.T = T;
            ut1_run.Beta = 0;
            ut1_run.H = Inf;

            ut1_run.WriteRun;
            ut1_run.Run;

            ut1_result = WamitResult(ut1_run);
            ut1_result.ReadResult;

            hydroComp = HydroBodyComp(ut1_result.HydroForces, cyl);

            RAO = hydroComp.Motions;
            figure;
            plot(T, abs(RAO));
            title({'WamitUT - testRAO', 'Heaving Cylinder'});
            ylabel('RAO');
            xlabel('T (s)');
            
            raoExp = [0, 0, 0, 0, 0, 0, 0, 0, 0, ...
                0.0001, 0.0004, 0.0014, ...
                0.0036, 0.0078, 0.0154, 0.0277, 0.0469, 0.0754, 0.1169, ...
                0.1764, 0.2617, 0.3849, 0.5653, 0.8333, 1.2251, 1.7252, ...
                2.1194, 2.1603, 1.9818, 1.7813, 1.6190, 1.4963, 1.4041, ...
                1.3337, 1.2792, 1.2361, 1.2017, 1.1738, 1.1508, 1.1317, ...
                1.1157, 1.1022, 1.0907, 1.0808, 1.0723, 1.0649, 1.0585, ...
                1.0528, 1.0479, 1.0435, 1.0396, 1.0362, 1.0331, 1.0303, ...
                1.0278, 1.0256, 1.0236, 1.0218, 1.0202, 1.0187]';
            
            testCase.verifyEqual(abs(RAO(9:end)), raoExp(9:end), 'AbsTol', 1e-2);
        end
        
        function testWF(testCase)
            name = 'ut2';
            folder = [WamitUT.getMainPath name];

            ut2_run = WamitRunCondition(folder, name);
            ut2_run.Rho = WamitUT.getRho;
            ut2_run.FieldArray = BemFieldArray([-20 -20 0], [0.4 0.4 1], [101 101 1]);

            cyl = WamitUT.createCylFB();
            cyl.Handle = 'Single';

            ut2_run.FloatingBodies = cyl;
            ut2_run.T = 2.8;
            ut2_run.Beta = 0;
            ut2_run.H = Inf;

            ut2_run.WriteRun;
            ut2_run.Run;

            ut2_result = WamitResult(ut2_run);
            ut2_result.ReadResult;

            wave = ut2_result.WaveArray;
            [X, Y] = wave.FieldPoints;

            % Radiated 1
            etaR1 = wave.Elevation('Radiated');

            figure; 
            subplot(2,2,1);
            pcolor(X,Y,real(etaR1{1}));
            shading flat; axis equal; axis tight; colorbar;
            set(gca, 'clim', [-0.03 0.03]);
            title({'WamitUT - testWF (T = 2.8)','','Radiated Wave - Real, Uncoupled'});

            % Radiated 2
            hydroComp = HydroBodyComp(ut2_result.HydroForces, cyl);

            wave.BodyMotions = hydroComp.Motions;

            etaR2 = wave.Elevation('Radiated');

            subplot(2,2,2);
            pcolor(X,Y,real(etaR2{1}));
            shading flat; axis equal; axis tight; colorbar;
            set(gca, 'clim', [-0.03 0.03]);
            title('Radiated Wave - Real, Coupled');

            % Diffracted
            etaD = wave.Elevation('Diffracted');

            subplot(2,2,3); 
            pcolor(X,Y,abs(etaD{1}));
            shading flat; axis equal; axis tight; colorbar;
            title('Diffracted Wave - Absolute Value');
            set(gca,'clim',[0.9 1.1]);

            % Total
            etaT = wave.Elevation('Total');

            subplot(2,2,4); 
            pcolor(X,Y,abs(etaT{1}));
            shading flat; axis equal; axis tight; colorbar;
            title('Total Wave - Absolute Value');
            set(gca,'clim',[0.9 1.1]);
        end
        
        function testSpectrum(testCase)
            Tp = 2.8;
            T = 5:-0.5:1;
            f = 1./T;
            Hs = 1;
            S = bretschneider(Hs, 1./Tp, f);
            Spec = WaveSpectrum(S,f);

            name = 'ut3';
            folder = [WamitUT.getMainPath name];

            ut3_run = WamitRunCondition(folder, name);
            ut3_run.Rho = WamitUT.getRho;
            ut3_run.FieldArray = BemFieldArray([-20 -20 0], [0.4 0.4 1], [101 101 1]);

            cyl = WamitUT.createCylFB();
            cyl.Handle = 'Spectrum';

            ut3_run.FloatingBodies = cyl;

            ut3_run.T = T;
            ut3_run.Beta = 0;
            ut3_run.H = Inf;

            ut3_run.WriteRun;
            ut3_run.Run;

            ut3_result = WamitResult(ut3_run);
            ut3_result.ReadResult;

            wave = ut3_result.WaveArray;
            [X, Y] = wave.FieldPoints;

            wave.IncidentWaveAmps = Spec.Amplitudes;

            hydroComp = HydroBodyComp(ut3_result.HydroForces, cyl);

            wave.BodyMotions = hydroComp.Motions;

            HsR = wave.SigWaveHeight('Radiated');

            figure; 
            subplot(3,1,1);
            pcolor(X,Y,HsR);
            shading flat; axis equal; axis tight; colorbar;
            title({'WamitUT - testSpectrum','','Radiated Wave - Hs'});

            HsD = wave.SigWaveHeight('Diffracted');

            subplot(3,1,2); 
            pcolor(X,Y,HsD);
            shading flat; axis equal; axis tight; colorbar;
            title('Diffracted Wave - Hs');
            set(gca, 'clim', [0.9 1.1])

            HsT = wave.SigWaveHeight('Total');

            subplot(3,1,3); 
            pcolor(X,Y,HsT);
            shading flat; axis equal; axis tight; colorbar;
            title('Total Wave - Hs');
            set(gca, 'clim', [0.9 1.1])
        end
        
        function test2cyl(testCase)
            tic;
            name = 'ut4';
            folder = [WamitUT.getMainPath name];

            ut4_run = WamitRunCondition(folder, name);
            ut4_run.Rho = WamitUT.getRho;
            ut4_run.FieldArray = BemFieldArray([-20 -20 0], [0.4 0.4 1], [101 101 1]);

            cyl1 = WamitUT.createCylFB();
            cyl2 = FloatingBody(cyl1);
            
            cyl1.Handle = 'Forward';
            cyl1.XYpos = [-2.5 0];

            cyl2.Handle = 'Aft';
            cyl2.XYpos = [2.5 0];

            ut4_run.FloatingBodies = [cyl1, cyl2];

            ut4_run.T = 2.8;
            ut4_run.Beta = 0;
            ut4_run.H = Inf;

            ut4_run.WriteRun;
            ut4_run.Run;
            toc

            tic
            ut4_result = WamitResult(ut4_run);
            ut4_result.ReadResult;
            toc

            wave = ut4_result.WaveArray;
            [X, Y] = wave.FieldPoints;

            hydroComp = HydroBodyComp(ut4_result.HydroForces, [cyl1 cyl2]);

            wave.BodyMotions = hydroComp.Motions;

            etaR = wave.Elevation('Radiated');

            figure; 
            subplot(2,2,1);
            pcolor(X,Y,abs(etaR{1}));
            shading flat; axis equal; axis tight; colorbar;
            set(gca, 'clim', [0 0.05]);
            title({'WamitUT - test2cyl, T = 2.8, Abs Elevation','','Forward Body, Radiated Wave'});

            subplot(2,2,2);
            pcolor(X,Y,abs(etaR{2}));
            shading flat; axis equal; axis tight; colorbar;
            set(gca, 'clim', [0 0.05]);
            title('Aft Body, Radiated Wave');

            etaD = wave.Elevation('Diffracted');

            subplot(2,2,3); 
            pcolor(X,Y,abs(etaD{1}));
            shading flat; axis equal; axis tight; colorbar;
            title('Diffracted Wave');
            set(gca, 'clim', [0.9 1.1]);

            etaT = wave.Elevation('Total');

            subplot(2,2,4); 
            pcolor(X,Y,abs(etaT{1}));
            shading flat; axis equal; axis tight; colorbar;
            title('Total Wave');
            set(gca, 'clim', [0.9 1.1]);
        end
        
        function testAttenRAO(testCase)
            name = 'ut6';
            folder = [WamitUT.getMainPath name];
            T = 0.1:0.1:6;
            rho = WamitUT.getRho;

            atten = WamitUT.createAttenFB();
            atten.Handle = 'RAO';

            ut6_run = WamitRunCondition(folder, name);
            ut6_run.Rho = rho;

            ut6_run.FloatingBodies = atten;
            ut6_run.T = T;
            ut6_run.Beta = 0;
            ut6_run.H = Inf;

            ut6_run.WriteRun;
            ut6_run.Run;

            ut6_result = WamitResult(ut6_run);
            ut6_result.ReadResult;
            
            hydroBody = HydroBodyComp(ut6_result.HydroForces, atten);
            d(2,2) = 3e4;
            hydroBody.SetDpto(d);
            P = hydroBody.Power;
            ef = IWaves.UnitEnergyFlux(rho, T, ut6_run.H);
            RAO = hydroBody.Motions;

            figure; 
            subplot(2,1,1);
            Ax = plotyy(T,abs(squeeze(RAO(:,1,1))),T,180/pi*abs(squeeze(RAO(:,1,2))));
            set(get(Ax(1),'Ylabel'),'String','RAO - heave') 
            set(get(Ax(2),'Ylabel'),'String','RAO - flex (deg)') 
            title({'WamitUT - testAttenRAO','','Attenuator Motions'});
            xlabel('T (s)');

            subplot(2,1,2);
            beam = 2*atten.Radius;
            plot(T, squeeze(P(:,1,2))./ef'./beam);
            title('Power');
            xlabel('T (s)');
            ylabel('RCW');
        end
        
        function testHingeModes(testCase)
            name = 'ut8';
            folder = [WamitUT.getMainPath name];

            rho = 1000;     
            
            len = 60;
            dia = 10;
            hingePos = [-10; 10];   
            sphereRad = 0.1*len;
            Nx = 80;              
            Ntheta = 24;           
            wec = FloatingSphereEndCyl(rho, len, dia/2, sphereRad, hingePos, Nx, Ntheta, 'Notch');  
            
            % Let's look at the Geometry..
            figure;
            plot(wec.PanelGeo);
            axis equal;
            title('WamitUT - testHingeModes Panelization');
            
            wec.Handle = 'atten';
            wec.Modes = ModesOfMotion([1 1 1 1 1 1 1 1]);   
            
            wec.WamIGenMds;
            
            wam_run = WamitRunCondition(folder, name);  

            wam_run.Rho = rho;      
            wam_run.T = 4:0.5:12;   
            wam_run.Beta = 0;       
            wam_run.H = Inf;        
            wam_run.FloatingBodies = wec;       

            wam_run.WriteRun;                   
            wam_run.Run;                           
            
            wam_result = WamitResult(wam_run);  
            wam_result.ReadResult;              
            hydroForces = wam_result.HydroForces;
            
            % Here, we have 8 DoF
            A = hydroForces.A;      
            B = hydroForces.B;      
            C = hydroForces.C;      
            Fex = hydroForces.Fex;  

            T = hydroForces.T;     
            figure;
            subplot(3,1,1);
            % DoF 7 and 8 are the hinge modes (or flex)
            % Because of symmetry the radiation forces 7-7 and 8-8 are the same
            axe = plotyy(T, [squeeze(A(:,7,7)) squeeze(A(:,8,8))], T, ...
                squeeze(A(:,3,7)));
            title({'WamitUT - testHingeModes', 'Added Mass'});
            ylabel(axe(1), 'kg*m^2');
            ylabel(axe(2), 'kg*m');
            legend('flex 1-flex 1', 'flex 2-flex 2', 'heave-flex 1');
            subplot(3,1,2);
            axe = plotyy(T, [squeeze(B(:,7,7)) squeeze(B(:,8,8))], T, ...
                squeeze(B(:,3,7)));
            title('Hydrodynamic Damping');
            ylabel(axe(1), 'Ns/m*m^2');
            ylabel(axe(2), 'Ns/m*m');
            legend('flex 1-flex 1', 'flex 2-flex 2', 'heave-flex 1');
            subplot(3,1,3);
            plot(T, abs([squeeze(Fex(:,1,7)) squeeze(Fex(:,1,8))]));
            title('Excitation Force');
            legend('flex 1', 'flex 2');
            xlabel('Period (s)');
            ylabel('Nm');
            
            hcomp = HydroBodyComp(hydroForces, wec);
            hcomp.C
        end
        
        function testReadFolder(testCase)
            name = 'ut7';
            folder = [WamitUT.getMainPath name];

            ut7_result = WamitResult;
            ut7_result.Folder = folder;
            ut7_result.RunName = name;
            ut7_result.ReadResult;
            
            cyl = WamitUT.createCylFB();

            hydroComp = HydroBodyComp(ut7_result.HydroForces, cyl);

            wave = ut7_result.WavePoints;
            points = wave.FieldPoints;
            wave.BodyMotions = hydroComp.Motions;

            % Total
            etaT = wave.Elevation('Total');

            figure; 
            scatter(points(:,1), points(:,2),4,abs(etaT{1}));
            shading flat; axis equal; axis tight; colorbar;
            title({'WamitUT - testReadFolder','','Total Wave - T = 2.8, Absolute Value'});
            set(gca,'clim',[0.9 1.1]);
        end
        
        function testZeroInfFreq(testCase)
            name = 'ut9';
            folder = [WamitUT.getMainPath name];
            T = [1 3 6 10];
           

            cyl = WamitUT.createCylFB();
            cyl.Modes = ModesOfMotion([1 0 1 0 1 0]);
            cyl.Handle = 'testInfZeroFreq';

            ut9_run = WamitRunCondition(folder, name);
            ut9_run.Rho = WamitUT.getRho;
            
            ut9_run.IncZeroFreq = true;
            ut9_run.IncInfFreq = true;

            ut9_run.FloatingBodies = cyl;
            ut9_run.T = T;
            ut9_run.Beta = 0;
            ut9_run.H = Inf;

            ut9_run.WriteRun;
            ut9_run.Run;

            ut1_result = WamitResult(ut9_run);
            ut1_result.ReadResult;

            hydroComp = HydroBodyComp(ut1_result.HydroForces, cyl);

            RAO = hydroComp.Motions;
            RAO = squeeze(RAO(:,1,3));
            figure;
            plot(T, abs(RAO));
            title({'WamitUT - testRAO', 'Heaving Cylinder'});
            ylabel('RAO');
            xlabel('T (s)');
        end
    end
    
    methods (Static, Access = private)
        
        function [cyl] = createCylFB()
            D = 500; 
            rho = WamitUT.getRho;

            diameter = 1;
            draft = 1.5;
            height = draft;
            
            Ntheta = 48;
            Nr = 8;
            Nz = 24;
            
            cyl = FloatingCylinder(rho, diameter/2, height, draft, Ntheta, Nr, Nz, 'NoInt');

            damping = zeros(6,6);
            damping(3,3) = D;
            cyl.Dpto = damping;
            cyl.Modes = ModesOfMotion([0 0 1 0 0 0]);
            cyl.XYpos = [0 0];
        end
        
        function [atten] = createAttenFB()
            rho = WamitUT.getRho;
            
            conePct = 0.05;
            len = 10;
            beam = sqrt(6/13);
            radius = beam/2;

            Nx = 80;
            Ntheta = 16;

            atten = FloatingAttenuator(rho, len, radius, conePct, Nx, Ntheta);  
            atten.Modes = ModesOfMotion([0 0 1 0 0 0 1]); 
        end
        
        function [mainPath] = getMainPath()
            mainPath = [mwavePath '\UnitTests\WamitRuns\'];
        end
        
        function [rho] = getRho()
            rho = 1000;
        end
        
    end
end