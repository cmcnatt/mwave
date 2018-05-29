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
classdef StlVideoUT <  matlab.unittest.TestCase
    
    methods (Test)
        function test1(testCase)            
            
            T = 7.3;
            h = Inf;
            beta = 0;
            lam = IWaves.T2Lam(T, h);
            slope = 50;
            a = lam/slope;
            
            dt = 0.1;
            time = 0:dt:(1*T-0.1);
            
            [cyl, mot] = StlVideoUT.getCylinderAndMotions(T, h, a, beta, time);
            
            stlCyl = cyl.PanelGeo.MakeStl;
            
            stlWave = StlWave(a, 2*pi/T, beta, h);
            stlWave.X = -50:2:50;
            stlWave.Y = stlWave.X;
            
            stlVideo = StlVideo;
            stlVideo.Bodies = stlCyl;
            stlVideo.Wave = stlWave;
            stlVideo.Time = time;
            stlVideo.Motions = mot;
            
            mov = movie(stlVideo, 'view', [10 10]);
            
            vidFile = [mwavePath '\UnitTests\files\stlVid_vid1.avi'];
            v = VideoWriter(vidFile);
            v.FrameRate = 1/dt;
            open(v);
            writeVideo(v, mov);
            close(v);
        end
        
        function test2(testCase)            
            
            T = 2:0.5:12;
            h = Inf;
            beta = 0;
            
            Hs = 2;
            Tp = 7.3;
            spec = Bretschneider(Hs, Tp, T);
            a = spec.Amplitudes('randPhase');
            
            dt = 0.1;
            time = 0:dt:(60-dt);
            
            [cyl, mot] = StlVideoUT.getCylinderAndMotions(T, h, a, beta, time);
            
            stlCyl = cyl.PanelGeo.MakeStl;
            
            stlWave = StlWave(a, 2*pi./T, beta, h);
            stlWave.X = -50:2:50;
            stlWave.Y = stlWave.X;
            
            stlVideo = StlVideo;
            stlVideo.Bodies = stlCyl;
            stlVideo.Wave = stlWave;
            stlVideo.Time = time;
            stlVideo.Motions = mot;
            
            mov = movie(stlVideo, 'view', [10 10]);
            
            vidFile = [mwavePath '\UnitTests\files\stlVid_vid2.avi'];
            v = VideoWriter(vidFile);
            v.FrameRate = 1/dt;
            open(v);
            writeVideo(v, mov);
            close(v);
        end
        
        function test3(testCase)            
            
            T = 2:0.5:12;
            h = Inf;
            beta = -pi:pi/45:pi;
            
            Hs = 2;
            Tp = 7.3;
            spec = Bretschneider(Hs, Tp, T, 'Cos2s', 2, 0, beta);
            a = spec.Amplitudes('randPhase');
            
            dt = 0.1;
            time = 0:dt:(60-dt);
            
            [cyl, mot] = StlVideoUT.getCylinderAndMotions(T, h, a, beta, time);
            
            stlCyl = cyl.PanelGeo.MakeStl;
            
            [Beta, Omega] = meshgrid(beta, 2*pi./T);
            
            stlWave = StlWave(a(:), Omega(:), Beta(:), h);
            stlWave.X = -50:2:50;
            stlWave.Y = stlWave.X;
            
            stlVideo = StlVideo;
            stlVideo.Bodies = stlCyl;
            stlVideo.Wave = stlWave;
            stlVideo.Time = time;
            stlVideo.Motions = mot;
            
            mov = movie(stlVideo, 'view', [10 10]);
            
            vidFile = [mwavePath '\UnitTests\files\stlVid_vid3.avi'];
            v = VideoWriter(vidFile);
            v.FrameRate = 1/dt;
            open(v);
            writeVideo(v, mov);
            close(v);
        end
    end
    
    methods (Static)
        function [wamRun] = getCylinderWamRun()
            
            folder = [mwavePath '\UnitTests\WamitRuns\tempRun'];  
            rho = 1000; 
            
            diameter = 10;
            draft = 10;
            height = 15;            
            
            Ntheta = 24;            
            Nr = 5;                 
            Nz = 12;                

            cyl = FloatingCylinder(rho, diameter/2, height, draft, Ntheta, Nr, Nz);
            cyl.Handle = 'Cylinder';                                           
            cyl.Modes = ModesOfMotion([1 1 1 1 1 1]);   
            
            wamRun = WamitRunCondition(folder, 'stlCyl');  

            wamRun.Rho = rho;         
            wamRun.Beta = 0;                 

            wamRun.FloatingBodies = cyl;       
            wamRun.WriteRun;                               
        end
        
        function [cyl, mot] = getCylinderAndMotions(T, h, a, beta, time)
            if nargin < 1
                T = 7.3;
            end
            if nargin < 2
                h = Inf;
            end
            if nargin < 3
                a = 1;
            end
            if nargin < 4
                beta = 0;
            end
            if nargin < 5
                time = 0:0.1:T(end)*5;
            end
            
            motFile = [mwavePath '\UnitTests\files\stlVid_mot0'];
            
            if T(1) == 7.3 && beta(1) == 0 && isinf(h) && exist([motFile '.mat'])
                load(motFile);
            else
                wamRun = StlVideoUT.getCylinderWamRun;
                cyl = wamRun.FloatingBodies;

                wamRun.T = T;
                wamRun.Beta = beta;
                wamRun.H = h;
                wamRun.WriteRun;   
                wamRun.Run;

                wamResult = WamitResult(wamRun);
                wamResult.ReadResult; 

                hydroForces = wamResult.FreqDomForces;
                hydroComp = FreqDomComp(hydroForces, cyl);

                Dpto = zeros(6,6);
                Dpto(3,3) = 100*10^3;
                hydroComp.SetDpto(Dpto);

                mot0 = hydroComp.Motions;
                
                if T(1) == 7.3 && beta(1) == 0 && isinf(h)
                    save(motFile, 'cyl', 'mot0');
                end
            end
            
            omega = 2*pi./T;
            Nt = length(T);
            Nb = length(beta);
            if any(size(a) ~= [Nt, Nb])
                error('The size of amplitudes must be the same as (number of periods) x (number of directions)');
            end
            
            if iscolumn(time)
                time = time';
            end
            
            mot = zeros(length(time), 6);
            for m = 1:Nt
                for n = 1:Nb
                    mot = mot + (real(a(m, n)*squeeze(mot0(m,n,:))*exp(1i*omega(m)*time))).';
                end
            end
        end
    end
end