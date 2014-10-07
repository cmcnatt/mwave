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
classdef MotionsVideo < handle
    
    properties (Access = private)
        t;
        h;
        omega;
        pointMot;
        wave;
        xWave;
        xlim; 
        zlim;
        motions;
        motionFuncs;
        bodyPoints;
        dof;
        showWave;
        isComp;
    end

    properties (Dependent)
        T;
        H;
        Motions;
        MotionFuncs;
        BodyPoints;
        DoF;
        ShowWave;
        Xlim;
        Zlim;
    end
    
    methods

        function [mv] = MotionsVideo(period, motions, motionFuncs, bodyPoints)
            
            if (length(period) ~= 1)
                error('Only a single period allowed');
            end
            
            if (period <= 0)
                error('Period must be a positive number');
            end
            
            mv.t = period;
            mv.omega = 2*pi/period;
            dof_ = length(motions);
            mv.dof = dof_;
            
            if (length(motionFuncs) ~= dof_)
                error('The number of motion functions must be the same as the number of motions (i.e. the DoF)');
            end
            
            mv.motions = motions;
            mv.motionFuncs = motionFuncs;
            
            [~, col] = size(bodyPoints);
            
            if (col ~= 3)
                error('The body points must be an Np x 3 vector of [x, y, z] body points');
            end
            
            mv.bodyPoints = bodyPoints;       
            
            mv.showWave = false;
        end
                
        function [t_] = get.T(mv)
            t_ = mv.t;
        end
        function [mv] = set.T(mv, t_)
            if (t_ <= 0)
                error('Period must be a positive value');
            end
            mv.isComp = false;
            mv.t = t_;
        end
        
        function [h_] = get.H(mv)
            h_ = mv.h;
        end
        function [mv] = set.H(mv, h_)
            if (h_ <= 0)
                error('The water depth must be positive');
            end
            mv.isComp = false;
            mv.h = h_;
        end
        
        function [mot] = get.Motions(mv)
            mot = mv.motions;
        end
        function [mv] = set.Motions(mv, mot)
            if (length(mot) ~= mv.dof)
                error('The number of motions must be the same as the DoF');
            end
            mv.isComp = false;
            mv.motions = mot;
        end
        
        function [mf] = get.MotionFuncs(mv)
            mf = mv.motionFuncs;
        end
        function [mv] = set.MotionFuncs(mv, motFuncs)
            if (length(motFuncs) ~= mv.dof)
                error('The number of motion functions must be the same as the DoF');
            end
            mv.isComp = false;
            mv.motionFuncs = motFuncs;
        end
        
        function [bp] = get.BodyPoints(mv)
            bp = mv.bodyPoints;
        end
        function [mv] = set.BodyPoints(mv, bodyPs)
            [~, col] = size(bodyPs);
            
            if (col ~= 3)
                error('The body points must be an Np x 3 vector of [x, y, z] body points');
            end
            mv.isComp = false;
            mv.bodyPoints = bodyPs;
        end
        
        function [do] = get.DoF(mv)
            do = mv.dof;
        end
        
        function [sw] = get.ShowWave(mv)
            sw = mv.showWave;
        end
        function [mv] = set.ShowWave(mv, sw)
            if (~isBool(sw))
                error('The ShowWave value must be boolean');
            end
            mv.isComp = false;
            mv.showWave = sw;
        end
        
        function [xli] = get.Xlim(mv)
            if (~mv.isComp)
                mv.computeMotWave;
            end
            xli = mv.xlim;
        end
        
        function [zli] = get.Zlim(mv)
            if (~mv.isComp)
                mv.computeMotWave;
            end
            zli = mv.zlim;
        end
        
        function [bodPos, wav, xwav] = GetPoints(mv, ti)
            
            if (~mv.isComp)
                mv.computeMotWave;
            end
            
            time = exp(1i*mv.omega*ti(1));
            
            if (mv.showWave)
                wav = real(mv.wave*time);
                xwav = mv.xWave;
            else
                wav = [];
                xwav = [];
            end
            
            bodPos = mv.bodyPoints + real(mv.pointMot*time);
        end
        
        function [mov] = movie(runTime, dt, motionsVids, varargin)
            % overloaded movie function
            
            [opts, args] = checkOptions({{'EqualZ'}, {'Title', 1}, {'YLabel', 1}}, varargin);
            
            equalz = opts(1);
            if (opts(2))
                titl = args{2}{1};
            else
                titl = -1;
            end
            
            if (opts(3))
                ylab = args{3}{1};
            else
                ylab = -1;
            end
                        
            ti = 0:dt:runTime;
            Nti = length(ti);
            
            Nmvs = length(motionsVids);
            
            figure;
            
            limx = zeros(Nmvs, 1);
            limz = zeros(Nmvs, 1);
            maxz = 0;
            for n = 1:Nmvs
                %plotaxes(n) = subplot(Nmvs, 1, n);
                
                axis equal;
                limx(n) = motionsVids(n).Xlim;
                limz(n) = motionsVids(n).Zlim;
                if (limz(n) > maxz)
                    maxz = limz(n);
                end
                %set(gca, 'xlim', [-limx limx], 'ylim', [-limz limz]);
            end
            
            if (equalz)
                limz = maxz*ones(Nmvs, 1);
            end
            
            nFrames = Nti;
            mov(1:nFrames) = struct('cdata',[], 'colormap',[]);
            
            for m = 1:Nti
                for n = 1:Nmvs
                    %set(gcf,'CurrentAxes',plotaxes(n));
                    subplot(Nmvs, 1, n);
                    cla;
                    set(gca, 'xlim', [-limx(n) limx(n)], 'ylim', [-limz(n) limz(n)]);
                    if (n == 1)
                        if (titl ~= -1)
                            title(titl);
                        end
                    end
                    
                    [bp, wav, xwav] = motionsVids(n).GetPoints(ti(m));
                    
                    if (motionsVids(n).ShowWave)
                        plot(xwav,wav);
                    end
                    hold on;
                    
                    if (iscell(ylab))
                        ylabel(ylab{n});
                    end
                    
                    scatter(bp(:,1), bp(:,3));
                    plot(bp(:,1), bp(:,3));
                    pause(0.5*dt);
                end
                mov(m) = getframe(gcf);
            end
        end
    end
    
    methods (Access = private)
        function [] = computeMotWave(mv)            
            mv.omega = 2*pi/mv.t;                      
                        
            Np = size(mv.bodyPoints, 1);
            
            mv.pointMot = zeros(Np, 3);
            
            for m = 1:Np
                for n = 1:mv.dof
                    mv.pointMot(m,:) = mv.pointMot(m,:) + mv.motions(n)*mv.motionFuncs(n).Evaluate(mv.bodyPoints(m,:));
                end
            end
            
            bp = mv.bodyPoints + abs(mv.pointMot);
            
            limx = max(abs(bp(:,1)));
            if (limx < 2)
                limx = 2;
            end
            
            mv.xlim = limx + 0.5*limx;
            
            limz = max(abs(bp(:,3)));
            if (limz < 1)
                limz = 1;
            end
            
            mv.zlim = limz + 0.1*limz;
            
            if (mv.showWave)
                if (isempty(mv.h))
                    error('Depth value (H) must be set to show the wave');
                end
                k = IWaves.SolveForK(mv.omega, mv.h, IWaves.G);
                mv.xWave = -mv.xlim:0.1:mv.xlim;
                mv.wave = exp(-1i*k*mv.xWave);
            end
            
            mv.isComp = true;
        end
    end
end