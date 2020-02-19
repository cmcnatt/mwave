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

classdef BodySurfWaveFieldUT < matlab.unittest.TestCase
    % Unit test for WaveCompClass

    methods (Test)
        
        function test0(testCase)
            % Just an RAO
            wam_run = WamitRunCondition(BodySurfWaveFieldUT.getFolder, 'test');  
            
            cyl = BodySurfWaveFieldUT.CreateCylFB;
            % cyl.GeoFile = 'cyl_mesh';  Commented out 21/01/2020 because
            % using programmatically create geometry
            % cyl.ISurfPan = 1;

            cyl.Modes = ModesOfMotion([1 0 1 0 1 0]);   

            wam_run.Rho = BodySurfWaveFieldUT.getRho;      
            wam_run.T = 4:0.2:10;   
            wam_run.Beta = 0;       
            wam_run.H = Inf;           

            wam_run.FloatingBodies = cyl;   

            wam_run.WriteRun;                   

            wam_run.Run;                          
            wam_result = WamitResult(wam_run);  
            wam_result.ReadResult;   

            comp = FreqDomComp(wam_result.FreqDomForces, cyl);
            xi = comp.Motions;
            
            figure;
            plot(comp.T, abs(xi(:,1,2)));
            title('BodySurfWaveFieldUT: test0: RAO');
            xlabel('T (s)');
            ylabel('Heave (m/m)');
        end
       
        function test1(testCase)
            % compute total wave field, and basic hydrostatic, use geo mesh
            wam_run = WamitRunCondition(BodySurfWaveFieldUT.getFolder, 'test');  
            
            cyl = BodySurfWaveFieldUT.CreateCylFB;
%             cyl.GeoFile = 'cyl_mesh_noInt';
%             cyl.ISurfPan = 0;

            cyl.Modes = ModesOfMotion([1 0 1 0 1 0]);   

            wam_run.Rho = BodySurfWaveFieldUT.getRho;      
            wam_run.T = 6;   
            wam_run.Beta = 0;       
            wam_run.H = Inf;           

            wam_run.FloatingBodies = cyl;   

            wam_run.ComputeBodyPoints = true;
            wam_run.ComputeVelocity = false;

            wam_run.WriteRun;                   

            wam_run.Run;                          
            wam_result = WamitResult(wam_run);  
            wam_result.ReadResult;   

            waveBody = wam_result.WaveBody;
            hydroForces = wam_result.FreqDomForces;
            comp = FreqDomComp(hydroForces, cyl);
            xi = comp.Motions;
            waveBody.BodyMotions = xi;

            
            press = waveBody.Pressure('Total');
            press = press{1};

            figure;
            subplot(2,1,1);
            surf(abs(press));
            axis equal
            title('BodySurfWaveFieldUT: test1: Panel Wave Body. Total Pressure');

            press = waveBody.Pressure('Hydrostatic');
            press = press{1};
            
            subplot(2,1,2);
            surf(press); 
            axis equal
            title('Hydrostatic Pressure');
        end
        
        function test2(testCase)
            % compute total wave field, and basic hydrostatic, use comp
            % mesh
            wam_run = WamitRunCondition(BodySurfWaveFieldUT.getFolder, 'test');  
            
            cyl = BodySurfWaveFieldUT.createCylFB;
            cyl.GeoFile = 'cyl_mesh';
            cyl.CompGeoFile = 'cyl_mesh_noInt';
            cyl.ISurfPan = 1;

            cyl.Modes = ModesOfMotion([1 0 1 0 1 0]);   

            wam_run.Rho = BodySurfWaveFieldUT.getRho;      
            wam_run.T = 6;   
            wam_run.Beta = 0;       
            wam_run.H = Inf;           

            wam_run.FloatingBodies = cyl;   

            wam_run.ComputeBodyPoints = true;
            wam_run.ComputeVelocity = false;

            wam_run.WriteRun;                   

            wam_run.Run;                          
            wam_result = WamitResult(wam_run);  
            wam_result.ReadResult;   

            waveBody = wam_result.WaveBody;
            hydroForces = wam_result.FreqDomForces;
            comp = FreqDomComp(hydroForces, cyl);
            xi = comp.Motions;
            waveBody.BodyMotions = xi;

            press = waveBody.Pressure('Total');
            press = press{1};

            figure;
            subplot(2,1,1);
            surf(abs(press));
            axis equal
            title('BodySurfWaveFieldUT: test2: Comp Wave Body. Total Pressure');

            press = waveBody.Pressure('Hydrostatic');
            press = press{1};
            
            subplot(2,1,2);
            surf(press); 
            axis equal
            title('Hydrostatic Pressure');
        end
        
        function test3(testCase)
            % compute total wave field, and basic hydrostatic, WAMIT hi order
            wam_run = WamitRunCondition(BodySurfWaveFieldUT.getFolder, 'test');  
            
            cyl = BodySurfWaveFieldUT.createCylFB;
            cyl.GeoFile = 'cyl_mesh_hi';
            cyl.WamILowHi = 1;
            cyl.WamPanelSize = 4;
            cyl.CompGeoFile = 'cyl_mesh_noInt';
            cyl.ISurfPan = 1;

            cyl.Modes = ModesOfMotion([1 0 1 0 1 0]);   

            wam_run.Rho = BodySurfWaveFieldUT.getRho;      
            wam_run.T = 6;   
            wam_run.Beta = 0;       
            wam_run.H = Inf;           

            wam_run.FloatingBodies = cyl;   

            wam_run.ComputeBodyPoints = true;
            wam_run.ComputeVelocity = false;

            wam_run.WriteRun;                   

            wam_run.Run;                          
            wam_result = WamitResult(wam_run);  
            wam_result.ReadResult;   

            waveBody = wam_result.WaveBody;
            hydroForces = wam_result.FreqDomForces;
            comp = FreqDomComp(hydroForces, cyl);
            xi = comp.Motions;
            waveBody.BodyMotions = xi;

            press = waveBody.Pressure('Total');
            press = press{1};

            figure;
            subplot(2,1,1);
            surf(abs(press));
            axis equal
            title('BodySurfWaveFieldUT: test3: Hi-order Wave Body. Total Pressure');

            press = waveBody.Pressure('Hydrostatic');
            press = press{1};
            
            subplot(2,1,2);
            surf(press); 
            axis equal
            title('Hydrostatic Pressure');
            
            testCase.assertFail('WAMIT Hi order method is not working currently');
        end
        
        function test4(testCase)
            % compute hydrostatic motions. Check calc
            wam_run = WamitRunCondition(BodySurfWaveFieldUT.getFolder, 'test');  
            
            cyl = BodySurfWaveFieldUT.createCylFB;
            cyl.GeoFile = 'cyl_mesh';
            cyl.CompGeoFile = 'cyl_mesh_noInt';
            cyl.ISurfPan = 1;

            cyl.Modes = ModesOfMotion([1 0 0 0 0 0]);   

            wam_run.Rho = BodySurfWaveFieldUT.getRho;      
            wam_run.T = 6;   
            wam_run.Beta = 0;       
            wam_run.H = Inf;           

            wam_run.FloatingBodies = cyl;   

            wam_run.ComputeBodyPoints = true;
            wam_run.ComputeVelocity = false;

            wam_run.WriteRun;                   

            wam_run.Run;                          
            wam_result = WamitResult(wam_run);  
            wam_result.ReadResult;   

            waveBody = wam_result.WaveBody;
            hydroForces = wam_result.FreqDomForces;
            comp = FreqDomComp(hydroForces, cyl);
            xi = comp.Motions;
            waveBody.BodyMotions = xi;

            press = waveBody.Pressure('MotionHydrostatic');
            press = press{1};
                        
            force = BodySurfWaveFieldUT.computeForceFromPress(waveBody.BodyGeo, waveBody.MotionFuncs, press);
            forceExp = hydroForces.C*xi;
            
            err = (force-forceExp)/forceExp*100;
            
            fprintf('\nForce \t\t\t= %4.2e,\nForce Expected \t= %4.2e\nError \t\t\t= %4.1f\n\n', force, forceExp, err)

            testCase.assertEqual(force, forceExp, 'RelTol', 0.01);
            
            figure;
            surf(abs(press));
            axis equal
            title('BodySurfWaveFieldUT: test4: MotionHydrostatic. surge');
        end
        
        function test5(testCase)
            % compute hydrostatic motions. Check calc - heave
            wam_run = WamitRunCondition(BodySurfWaveFieldUT.getFolder, 'test');  
            
            cyl = BodySurfWaveFieldUT.createCylFB;
            cyl.GeoFile = 'cyl_mesh';
            cyl.CompGeoFile = 'cyl_mesh_noInt';
            cyl.ISurfPan = 1;

            cyl.Modes = ModesOfMotion([0 0 1 0 0 0]);   

            wam_run.Rho = BodySurfWaveFieldUT.getRho;      
            wam_run.T = 6;   
            wam_run.Beta = 0;       
            wam_run.H = Inf;           

            wam_run.FloatingBodies = cyl;   

            wam_run.ComputeBodyPoints = true;
            wam_run.ComputeVelocity = false;

            wam_run.WriteRun;                   

            wam_run.Run;                          
            wam_result = WamitResult(wam_run);  
            wam_result.ReadResult;   

            waveBody = wam_result.WaveBody;
            hydroForces = wam_result.FreqDomForces;
            comp = FreqDomComp(hydroForces, cyl);
            xi = comp.Motions;
            waveBody.BodyMotions = xi;

            press = waveBody.Pressure('MotionHydrostatic');
            press = press{1};
            
            force = BodySurfWaveFieldUT.computeForceFromPress(...
                waveBody.BodyGeo, waveBody.MotionFuncs, press);
            forceExp = hydroForces.C*xi;
            
            err = (force-forceExp)/forceExp*100;
            
            fprintf('\nForce \t\t\t= %4.2e,\nForce Expected \t= %4.2e\nError \t\t\t= %4.1f\n\n', force, forceExp, err)

            testCase.assertEqual(force, forceExp, 'RelTol', 0.01)
            
            figure;
            surf(abs(press));
            axis equal
            title('BodySurfWaveFieldUT: test5: MotionHydrostatic. heave');
        end
        
        function test6(testCase)
            % compute hydrostatic motions. Check calc - pitch
            wam_run = WamitRunCondition(BodySurfWaveFieldUT.getFolder, 'test');  
            
            cyl = BodySurfWaveFieldUT.createCylFB;
            cyl.GeoFile = 'cyl_mesh0';
            cyl.CompGeoFile = 'cyl_mesh_noInt0';
            cyl.Cg = [0 0 0];
            cyl.ISurfPan = 1;

            cyl.Modes = ModesOfMotion([0 0 0 0 1 0]);   

            wam_run.Rho = BodySurfWaveFieldUT.getRho;      
            wam_run.T = 6;   
            wam_run.Beta = 0;       
            wam_run.H = Inf;           

            wam_run.FloatingBodies = cyl;   

            wam_run.ComputeBodyPoints = true;
            wam_run.ComputeVelocity = false;

            wam_run.WriteRun;                   

            wam_run.Run;                          
            wam_result = WamitResult(wam_run);  
            wam_result.ReadResult;   

            waveBody = wam_result.WaveBody;
            hydroForces = wam_result.FreqDomForces;
            comp = FreqDomComp(hydroForces, cyl);
            xi = comp.Motions;
            waveBody.BodyMotions = xi;

            press = waveBody.Pressure('MotionHydrostatic');
            press = press{1};
            
            force = BodySurfWaveFieldUT.computeForceFromPress(...
                waveBody.BodyGeo, waveBody.MotionFuncs, press);
            forceExp = hydroForces.C*xi;
            
            err = (force-forceExp)/forceExp*100;
            
            fprintf('\nForce \t\t\t= %4.2e,\nForce Expected \t= %4.2e\nError \t\t\t= %4.1f\n\n', force, forceExp, err)

            testCase.assertEqual(force, forceExp, 'RelTol', 0.01)
            
            figure;
            surf(abs(press));
            axis equal
            title('BodySurfWaveFieldUT: test6: MotionHydrostatic. pitch at 0');
        end
        
        function test7(testCase)
            % compute hydrostatic motions. Check calc - pitch
            wam_run = WamitRunCondition(BodySurfWaveFieldUT.getFolder, 'test');  
            
            cyl = BodySurfWaveFieldUT.createCylFB;
            cyl.GeoFile = 'cyl_mesh';
            cyl.CompGeoFile = 'cyl_mesh_noInt';
            cyl.ISurfPan = 1;

            cyl.Modes = ModesOfMotion([0 0 0 0 1 0]);   

            wam_run.Rho = BodySurfWaveFieldUT.getRho;      
            wam_run.T = 6;   
            wam_run.Beta = 0;       
            wam_run.H = Inf;           

            wam_run.FloatingBodies = cyl;   

            wam_run.ComputeBodyPoints = true;
            wam_run.ComputeVelocity = false;

            wam_run.WriteRun;                   

            wam_run.Run;                          
            wam_result = WamitResult(wam_run);  
            wam_result.ReadResult;   

            waveBody = wam_result.WaveBody;
            hydroForces = wam_result.FreqDomForces;
            comp = FreqDomComp(hydroForces, cyl);
            xi = comp.Motions;
            waveBody.BodyMotions = xi;

            press = waveBody.Pressure('MotionHydrostatic');
            press = press{1};
            
            force = BodySurfWaveFieldUT.computeForceFromPress(...
                waveBody.BodyGeo, waveBody.MotionFuncs, press);
            forceExp = hydroForces.C*xi;
            
            err = (force-forceExp)/forceExp*100;
            
            fprintf('\nForce \t\t\t= %4.2e,\nForce Expected \t= %4.2e\nError \t\t\t= %4.1f\n\n', force, forceExp, err)

            testCase.assertEqual(force, forceExp, 'RelTol', 0.02)
            
            figure;
            surf(abs(press));
            axis equal
            title('BodySurfWaveFieldUT: test7: MotionHydrostatic. pitch');
        end
        
        function test8(testCase)
            % compute hydrostatic motions. Check calc - total forces
            wam_run = WamitRunCondition(BodySurfWaveFieldUT.getFolder, 'test');  
            
            cyl = BodySurfWaveFieldUT.createCylFB;
            cyl.GeoFile = 'cyl_mesh';
            cyl.CompGeoFile = 'cyl_mesh_noInt';
            cyl.ISurfPan = 1;

            cyl.Modes = ModesOfMotion([1 0 1 0 1 0]);   

            wam_run.Rho = BodySurfWaveFieldUT.getRho;      
            wam_run.T = 6;   
            wam_run.Beta = 0;       
            wam_run.H = Inf;           

            wam_run.FloatingBodies = cyl;   

            wam_run.ComputeBodyPoints = true;
            wam_run.ComputeVelocity = false;

            wam_run.WriteRun;                   

            wam_run.Run;                          
            wam_result = WamitResult(wam_run);  
            wam_result.ReadResult;   

            waveBody = wam_result.WaveBody;
            hydroForces = wam_result.FreqDomForces;
            comp = FreqDomComp(hydroForces, cyl);
            xi = comp.Motions;
            waveBody.BodyMotions = xi;

            press = waveBody.Pressure('Hydrostatic');
            pressHs0 = press{1};
            
            press = waveBody.Pressure('MotionHydrostatic');
            pressHsM = press{1};
            
            press = waveBody.Pressure('Total');
            pressHd = press{1};
            
            geo = waveBody.BodyGeo;
            funcs = waveBody.MotionFuncs;
            
            forceHs0 = zeros(3, 1);
            forceHsM = zeros(3, 1);
            forceHd = zeros(3, 1);
            
            for n = 1:length(funcs)
                forceHs0(n) = BodySurfWaveFieldUT.computeForceFromPress(geo, funcs(n), pressHs0);
                forceHsM(n) = BodySurfWaveFieldUT.computeForceFromPress(geo, funcs(n), pressHsM);
                forceHd(n) = BodySurfWaveFieldUT.computeForceFromPress(geo, funcs(n), pressHd);    
            end
            
            xi = squeeze(xi);
            m = cyl.M(1,1);
            g = IWaves.G;
            forceHs0e = m*g*[xi(3); -1; 0];
            
            forceHsMe = hydroForces.C*xi;
            
            Fex = squeeze(hydroForces.Fex);
            A = squeeze(hydroForces.A);
            B = squeeze(hydroForces.B);
            w = 2*pi/hydroForces.T;
            forceHde = -(Fex + w^2*A*xi - 1i*w*B*xi);
            
            force = forceHs0 + forceHsM + forceHd;
            forceExp = forceHs0e + forceHsMe + forceHde;
            
            err = abs(force-forceExp)./abs(forceExp)*100;
            
            for n = 1:3
                fprintf('\nForce \t\t\t= %4.2e,\nForce Expected \t= %4.2e\nError \t\t\t= %4.1f\n\n', abs(force(n)), abs(forceExp(n)), err(n))
            end

            for n = 1:3
                testCase.assertEqual(force(n), forceExp(n), 'RelTol', 0.02)
            end
        end
        
        function test9(testCase)
            % compute hydrostatic motions. Check calc with built in test
            wam_run = WamitRunCondition(BodySurfWaveFieldUT.getFolder, 'test');  
            
            cyl = BodySurfWaveFieldUT.createCylFB;
            cyl.GeoFile = 'cyl_mesh';
            cyl.CompGeoFile = 'cyl_mesh_noInt';
            cyl.ISurfPan = 1;

            cyl.Modes = ModesOfMotion([1 0 1 0 1 0]);   

            wam_run.Rho = BodySurfWaveFieldUT.getRho;      
            wam_run.T = 6;   
            wam_run.Beta = 0;       
            wam_run.H = Inf;           

            wam_run.FloatingBodies = cyl;   

            wam_run.ComputeBodyPoints = true;
            wam_run.ComputeVelocity = false;

            wam_run.WriteRun;                   

            wam_run.Run;                          
            wam_result = WamitResult(wam_run);  
            wam_result.ReadResult;   

            waveBody = wam_result.WaveBody;
            hydroForces = wam_result.FreqDomForces;
            comp = FreqDomComp(hydroForces, cyl);
            xi = comp.Motions;
            waveBody.BodyMotions = xi;
            
            dofSur = 1;
            dofHea = 2;
            dofPit = 3;
            [forceHs0, forceHsM, forceHd, forceHs0e, forceHsMe, forceHde] ...
                = waveBody.CheckPressureForce(cyl.M(1,1), dofSur, dofHea, dofPit, hydroForces);
            
            force = forceHs0 + forceHsM + forceHd;
            forceExp = forceHs0e + forceHsMe + forceHde;
            
            for n = 1:3
                testCase.assertEqual(force(n), forceExp(n), 'RelTol', 0.02)
            end
        end
        
        function test10(testCase)
            % use field points as body points
            wam_run = WamitRunCondition(BodySurfWaveFieldUT.getFolder, 'test');  
            
            cyl = BodySurfWaveFieldUT.createCylFB;
            cyl.GeoFile = 'cyl_mesh';
            cyl.CompGeoFile = 'cyl_mesh_noInt';
            cyl.ISurfPan = 1;

            cyl.Modes = ModesOfMotion([1 0 1 0 1 0]);   

            wam_run.Rho = BodySurfWaveFieldUT.getRho;      
            wam_run.T = 6;   
            wam_run.Beta = 0;       
            wam_run.H = Inf;           

            wam_run.FloatingBodies = cyl;   

            wam_run.FieldPointAsBodyDistance = 0.05;
            wam_run.ComputeBodyPoints = true;
            wam_run.ComputeVelocity = false;

            wam_run.WriteRun;                   

            wam_run.Run;                          
            wam_result = WamitResult(wam_run);  
            wam_result.ReadResult;   

            waveBody = wam_result.WaveBody;
            hydroForces = wam_result.FreqDomForces;
            comp = FreqDomComp(hydroForces, cyl);
            xi = comp.Motions;
            waveBody.BodyMotions = xi;
            
            dofSur = 1;
            dofHea = 2;
            dofPit = 3;
            [forceHs0, forceHsM, forceHd, forceHs0e, forceHsMe, forceHde] ...
                = waveBody.CheckPressureForce(cyl.M(1,1), dofSur, dofHea, dofPit, hydroForces);
            
            force = forceHs0 + forceHsM + forceHd;
            forceExp = forceHs0e + forceHsMe + forceHde;
            
            for n = 1:3
                testCase.assertEqual(force(n), forceExp(n), 'RelTol', 0.02)
            end
        end
        
        function test11(testCase)
            % use field points as body points, and full cylinder
            wam_run = WamitRunCondition(BodySurfWaveFieldUT.getFolder, 'test');  
            
            cyl = BodySurfWaveFieldUT.createCylFB;
            cyl.GeoFile = 'cyl_mesh';
            cyl.CompGeoFile = 'cyl_mesh_F';
            cyl.ISurfPan = 1;

            cyl.Modes = ModesOfMotion([1 0 1 0 1 0]);   

            wam_run.Rho = BodySurfWaveFieldUT.getRho;      
            wam_run.T = 6;   
            wam_run.Beta = 0;       
            wam_run.H = Inf;           

            wam_run.FloatingBodies = cyl;   

            wam_run.FieldPointAsBodyDistance = 0.05;
            wam_run.ComputeBodyPoints = true;
            wam_run.ComputeVelocity = false;

            wam_run.WriteRun;                   

            wam_run.Run;                          
            wam_result = WamitResult(wam_run);  
            wam_result.ReadResult;   

            waveBody = wam_result.WaveBody;
            hydroForces = wam_result.FreqDomForces;
            comp = FreqDomComp(hydroForces, cyl);
            xi = comp.Motions;
            waveBody.BodyMotions = xi;
            
            dofSur = 1;
            dofHea = 2;
            dofPit = 3;
            [forceHs0, forceHsM, forceHd, forceHs0e, forceHsMe, forceHde] ...
                = waveBody.CheckPressureForce(cyl.M(1,1), dofSur, dofHea, dofPit, hydroForces);
            
            force = forceHs0 + forceHsM + forceHd;
            forceExp = forceHs0e + forceHsMe + forceHde;
            
            for n = 1:3
                testCase.assertEqual(force(n), forceExp(n), 'RelTol', 0.02)
            end
        end
        
        function test12(testCase)
            % write to file
            wam_run = WamitRunCondition(BodySurfWaveFieldUT.getFolder, 'test');  
            
            cyl = BodySurfWaveFieldUT.createCylFB;
            cyl.GeoFile = 'cyl_mesh';
            cyl.CompGeoFile = 'cyl_mesh_noInt';
            cyl.ISurfPan = 1;

            cyl.Modes = ModesOfMotion([1 0 1 0 1 0]);   

            wam_run.Rho = BodySurfWaveFieldUT.getRho;      
            wam_run.T = 6;   
            wam_run.Beta = 0;       
            wam_run.H = Inf;           

            wam_run.FloatingBodies = cyl;   

            wam_run.ComputeBodyPoints = true;
            wam_run.ComputeVelocity = false;

            wam_run.WriteRun;                   

            wam_run.Run;                          
            wam_result = WamitResult(wam_run);  
            wam_result.ReadResult;   

            waveBody = wam_result.WaveBody;
            hydroForces = wam_result.FreqDomForces;
            comp = FreqDomComp(hydroForces, cyl);
            xi = comp.Motions;
            waveBody.BodyMotions = xi;
            
            folder = [mwavePath '\UnitTests\writeFiles'];
            name = ['cylT' num2str(wam_run.T)];
            time = 0;
            A = 0.5;
            waveBody.WritePressureFile(folder, name, time, A)
        end
        
        function energyFlux(testCase)
            % compute energy flux body surf
            wam_run = WamitRunCondition(BodySurfWaveFieldUT.getFolder, 'test');  
            
            cyl = BodySurfWaveFieldUT.CreateCylFB;
%             cyl.GeoFile = 'cyl_mesh_noInt';
%             cyl.ISurfPan = 0;

            cyl.Modes = ModesOfMotion([1 0 1 0 1 0]);   

            wam_run.Rho = BodySurfWaveFieldUT.getRho;      
            wam_run.T = 6;   
            wam_run.Beta = 0;       
            wam_run.H = Inf;    

            wam_run.FloatingBodies = cyl;   

            wam_run.ComputeBodyPoints = true;
            wam_run.ComputeVelocity = true;

            wam_run.WriteRun;                   

            wam_run.Run;                          
            wam_result = WamitResult(wam_run);  
            wam_result.ReadResult;   

            waveBody = wam_result.WaveBody;
            hydroForces = wam_result.FreqDomForces;
            
            comp = FreqDomComp(hydroForces, cyl);

            waveBody.BodyMotions = comp.Motions;

            flux0 = waveBody.EnergyFlux([], 'Total');
            fluxD = waveBody.EnergyFlux([], 'Diffracted');
            flux0 = flux0{1};
            figure; 
            subplot(2,1,1);
            surf(flux0); 
            axis equal
            title('No PTO Damping');
            set(gca, 'view', [0 -40]);
            
            Dpto = zeros(3,3);
            Dpto(2,2) = 10^6;
            comp.SetDpto(Dpto);
            waveBody.BodyMotions = comp.Motions;
            
            fluxP = waveBody.EnergyFlux([], 'Total');
            fluxP = fluxP{1};
            subplot(2,1,2);
            surf(fluxP); 
            axis equal
            title('With PTO Damping');
            set(gca, 'view', [0 -40]);
            
            powComp = comp.Power;
            powComp = powComp(1,1,2);
            powSurf = -sum(fluxP.Values.*fluxP.Areas);
            
            testCase.assertEqual(powSurf, powComp, 'RelTol', 0.01);
        end
    end
    
    methods (Static, Access = public)
        
        function [cyl] = CreateCylFB()
            diameter = 10;
            draft = 10;
            height = 12; 
            
            Ntheta = 24;
            Nr = 8;
            Nz = 12;
            
            rho = BodySurfWaveFieldUT.getRho;
            
            cyl = FloatingCylinder(rho, diameter/2, height, draft, Ntheta, Nr, Nz);

        end
    end
    
    methods (Static, Access = private)
        
        function [cyl] = createCylFBOld()
            diameter = 10;
            draft = 10;
            height = 12;           

            cg = [0 0 -draft+height/2];           
            M = zeros(6,6);       
            
            rho = BodySurfWaveFieldUT.getRho;

            mass = rho*pi*(diameter/2)^2*draft; 
            Ixx = 1/12*mass*(3*(diameter/2)^2 + height^2);
            Iyy = Ixx;
            Izz = mass*(diameter/2)^2/2;

            M(1,1) = mass;
            M(2,2) = mass;
            M(3,3) = mass;
            M(4,4) = Ixx;
            M(5,5) = Iyy;
            M(6,6) = Izz;

            cyl = FloatingBody();           
            cyl.Handle = 'Cylinder';        
            cyl.M = M;                      
            cyl.Cg = cg;    
        end
        
        function [mainPath] = getFolder()
            mainPath = [mwavePath '\UnitTests\WamitRuns\tempRun'];
        end
        
        function [rho] = getRho()
            rho = 1025;
        end
        
        function [force] = computeForceFromPress(geo, func, press)
            cents = geo.Centroids;
            norms = geo.Normals;
            areas = geo.Areas;
            
            normMot = zeros(geo.Count, 3);
            dotNorm = zeros(geo.Count, 1);
            
            for n = 1:geo.Count
                normMot(n, :) = func.Evaluate(cents(n,:));
                dotNorm(n) = dot(normMot(n,:), norms(n,:));
            end

            force = sum(dotNorm.*areas.*press.Values);
        end
    end
end