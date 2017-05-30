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
classdef IEnergyCompUT < matlab.unittest.TestCase
    
    methods (Test)
        function testPowerRAO(testCase)
            
            comp = IEnergyCompUT.HingeBargeComp;
            pow = comp.PowerRAO;
            
            rho = 1025;
            wid = 4;
            CW = comp.PowerRAO('CW', rho);
            
            figure; 
            subplot(2,2,1:2)
            for n = 1:2
                geo = comp.Bodies(n).PanelGeo;
                geo.Translate([comp.Bodies(n).XYpos 0]);
                plot(geo);
                hold on;
            end
            axis equal
            
            title('IEnergyCompUT - testPowerRAO');
            
            subplot(2,2,3);
            plot(comp.T, pow);
            xlabel('T (s)');
            ylabel('Power RAO (kW/m^2)')
            
            subplot(2,2,4);
            plot(comp.T, CW./wid);
            xlabel('T (s)');
            ylabel('Relative Capture Width (m)')
        end
        
        function testAveragePower(testCase)
            rho = 1025;
            Hs = 2;
            Tp = 8;
            comp = IEnergyCompUT.HingeBargeComp;
            spec = Bretschneider(Hs, Tp, comp.T);
            
            avPow = comp.AveragePower(spec);
            
            avCW = comp.AveragePower(spec, 'CW', rho);
            
            fprintf('\n~~~~~~~~~~~~~\n\nBretschneider: Hs = %4.1f m, Tp = %4.1f s, Average Power = %4.1f kW\n\n~~~~~~~~~~~~\n', Hs, Tp, avPow);
            
            Ef = spec.EnergyFlux(rho, comp.H)./1000;  % energy flux in kW
            avCW2 = avPow/sum(Ef);
            
            testCase.verifyEqual(avCW, avCW2, 'AbsTol', 1e-9);
        end
        
        function testAnnualEnergyProd(testCase)
            
            comp = IEnergyCompUT.HingeBargeComp;
            
            waveClim = CreateClimates.EMEC(comp.T);
            AEP = comp.AnnualEnergyProd(waveClim);
            
            fprintf('\n~~~~~~~~~~~~~\n\nAEP at EMEC: %4.0f MWh\n\n~~~~~~~~~~~~\n', AEP);
            
            Nsamp = 100;
            AEPs = comp.AnnualEnergyProd(waveClim, 'Sample', [Nsamp 1]);
            
            testCase.verifyEqual(length(AEPs), Nsamp);
        end
    end
    
    methods(Static)
        function [comp] = HingeBargeComp()
            rho = 1000;
            lam = 6:2:120;
            h = Inf;
            T = IWaves.Lam2T(lam, h);
            
            lens = [15 15];
            wid = 4;
            hei = 3;
            draft = 2;           
            
            name = 'hb';
            folder = [mwavePath '\UnitTests\WamitRuns\tempRun'];
            mods = ModesOfMotion([1 1 1 1 1 1]); 

            Nx = 15;
            Ny = 4;
            Nz = 2;
            for n = 1:2
                bodies(n) = FloatingTrap(rho, lens(n), lens(n), 0, wid, ...
                    hei, draft, Nx, Ny, Nz);
                
                bodies(n).Modes = mods;
            end
            
            space = mean(lens)/20;
            bodies(1).XYpos = [-(lens(1) + space)/2 0];
            bodies(2).XYpos = [(lens(2) + space)/2 0];
            
            wam_run = WamitRunCondition(folder, name);  

            wam_run.Rho = rho;      
            wam_run.T = T;

            wam_run.Beta = 0; 
            wam_run.H = h;    
            
            wam_run.FloatingBodies = bodies;
            
            wam_run.MaxItt = 200;
            wam_run.WriteRun;            
            wam_run.Run;   

            wam_result = WamitResult(wam_run); 
            wam_result.ReadResult;
            
            forces = wam_result.FreqDomForces;
            
            cg1 = bodies(1).CgGlobal;
            cg2 = bodies(2).CgGlobal;

            hin = [0, 0, hei/2-draft];

            P = ConstraintMatComp.HingedBodies([cg1; cg2], hin);

            comp = FreqDomComp(forces, bodies, 'Constrained', P);
            
            Dpto = zeros(7,7);
            Dpto(7,7) = 1.8*10^6;
            comp.SetDpto(Dpto);
        end
    end
end