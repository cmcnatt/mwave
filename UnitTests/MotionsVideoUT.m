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
classdef MotionsVideoUT < matlab.unittest.TestCase

    methods (Test)
       
        function test1(testCase)
            load mv_ut1_comp;
            
            D = zeros(7,7);
            D(7,7) = 8e4;
            hydroComp.SetDpto(D);
           
            mot = hydroComp.Motions;
            mot = squeeze(mot);
            
            len = wec.Length;
            rad = wec.Radius;
            hingePos = wec.HingePos;
            
            T = hydroComp.T;
            h = hydroComp.H;
            
            funcs = wec.Modes.MotionFuncs;
            
            bodyPnts = [-len/2 0 0; hingePos 0 0; len/2 0 0; len/2 0 -rad; hingePos 0 -rad; -len/2 0 -rad; -len/2 0 0];
            
            mv = MotionsVideo(T, mot, funcs, bodyPnts);
            mv.ShowWave = true;
            mv.H = h;
            
            dt = 0.1;
            runT = 5*T;
            
            mov = movie(runT, dt, mv);
        end
    end
    
    methods (Static, Access = public)
        
        function [] = computeHydroMotions()
            
            T = 2.5;
            h = 10;
            len = 10;
            rho = 1000;
            dia = 1;
            hingePos = 1;
            sphereRad = 0.1*len;
            
            Nx = 80;
            Ntheta = 24;
            
            wec = FloatingSphereEndCyl(rho, len, dia/2, sphereRad, hingePos, Nx, Ntheta);  
            wec.Handle = 'atten';            
            
            name = 'mv_ut1';
            folder = [mwavePath 'UnitTests\WamitRuns\'];
            
            comp_run = WamitRunCondition(folder, name);
            comp_run.Rho = rho;
            comp_run.FloatingBodies = wec;
            comp_run.T = T;
            comp_run.Beta = 0;
            comp_run.H = h;
            
            comp_run.WriteRun;
            comp_run.RunWamit;
            
            comp_result = WamitResult(comp_run);
            comp_result.ReadResult;
            
            hydroComp = HydroBodyComp(comp_result.HydroForces, wec);
            
            save([mwavePath 'UnitTests\files\mv_ut1_comp'], 'wec', 'hydroComp');
        end
        
    end
end