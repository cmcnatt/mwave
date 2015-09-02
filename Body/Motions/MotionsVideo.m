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
        nBod;
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
        waveA;
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
        WaveAmp;
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
            
            if (iscell(motions))
                mv.nBod = length(motions);
                                
                if (~iscell(motionFuncs))
                    error('If motions are provided as a cell, then so must the motionFuncs and the bodyPoints');
                end
                
                if (~iscell(bodyPoints))
                    error('If motions are provided as a cell, then so must the motionFuncs and the bodyPoints');
                end
                
                if (mv.nBod == 1)
                    mv.motions = motions{1};
                    mv.motionFuncs = motionFuncs{1};
                    mv.bodyPoints = bodyPoints{1};  
                else
                    mv.motions = motions;
                    mv.motionFuncs = motionFuncs;
                    mv.bodyPoints = bodyPoints;  
                end
            else
                mv.nBod = 1;
                mv.motions = motions;
                mv.motionFuncs = motionFuncs;
                mv.bodyPoints = bodyPoints;  
            end
            
            dofs = zeros(mv.nBod,1);
            for n = 1:mv.nBod
                if (iscell(motions))
                    motn = motions{n};
                    motFunn = motionFuncs{n};
                    bodPntn = bodyPoints{n};
                else
                    motn = motions;
                    motFunn = motionFuncs;
                    bodPntn = bodyPoints;
                end
                
                dofs(n) = length(motn);
                
                if (length(motFunn) ~= dofs(n))
                    error('The number of motion functions must be the same as the number of motions (i.e. the DoF)');
                end
                
                [~, col] = size(bodPntn);
                
                if (col ~= 3)
                    error('The body points must be an Np x 3 vector of [x, y, z] body points');
                end
            end
            
            mv.dof = dofs;
            mv.waveA = 1;  
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
        
        function [a_] = get.WaveAmp(mv)
            a_ = mv.waveA;
        end
        function [mv] = set.WaveAmp(mv, a_)
            if (a_ <= 0)
                error('The wave amplitude must be positive');
            end
            mv.isComp = false;
            mv.waveA = a_;
        end
        
        function [mot] = get.Motions(mv)
            mot = mv.motions;
        end
        function [mv] = set.Motions(mv, mot)
            if (mv.nBod > 1)
                if (~iscell(mot) || (length(mot) ~= mv.nBod))
                    error('For multiple bodies, input must be cell of length equal to the number of bodies');
                end
            end
            for n = 1:mv.nBod
                if (iscell(mot))
                    mo = mot{n};
                else
                    mo = mot;
                end
                if (length(mo) ~= mv.dof(n))
                    error('The number of motions must be the same as the DoF');
                end
            end
            mv.isComp = false;
            mv.motions = mot;
        end
        
        function [mf] = get.MotionFuncs(mv)
            mf = mv.motionFuncs;
        end
        function [mv] = set.MotionFuncs(mv, motFuncs)
            if (mv.nBod > 1)
                if (~iscell(motFuncs) || (length(motFuncs) ~= mv.nBod))
                    error('For multiple bodies, input must be cell of length equal to the number of bodies');
                end
            end
            for n = 1:mv.nBod
                if (iscell(mot))
                    mf = motFuncs{n};
                else
                    mf = motFuncs;
                end
                if (length(mf) ~= mv.dof(n))
                    error('The number of motion functions must be the same as the DoF');
                end
            end
            mv.isComp = false;
            mv.motionFuncs = motFuncs;
        end
        
        function [bp] = get.BodyPoints(mv)
            bp = mv.bodyPoints;
        end
        function [mv] = set.BodyPoints(mv, bodyPs)
            if (mv.nBod > 1)
                if (~iscell(bodyPs) || (length(bodyPs) ~= mv.nBod))
                    error('For multiple bodies, input must be cell of length equal to the number of bodies');
                end
            end
            
            for n = 1:mv.nBod
                if (iscell(bodyPs))
                    bps = bodyPs{n};
                else
                    bps = bodyPs;
                end
                [~, col] = size(bps);

                if (col ~= 3)
                    error('The body points must be an Np x 3 vector of [x, y, z] body points');
                end
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
        function [mv] = set.Xlim(mv, xli)
            if (length(xli) ~= 1)
                error('Xlim must be a scalar. Limits are +/-Xlim');
            end
            mv.xlim = xli;
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
            
            if (mv.nBod == 1)
                bodPos = mv.bodyPoints + real(mv.pointMot*time);
            else
                bodPos = cell(mv.nBod, 1);
                for n = 1:mv.nBod
                    bodPos{n} = mv.bodyPoints{n} + real(mv.pointMot{n}*time);
                end
            end
        end
        
        function [mov] = movie(runTime, dt, motionsVids, varargin)
            % overloaded movie function
            
            [opts, args] = checkOptions({{'EqualZ'}, {'Title', 1}, ...
                {'YLabel', 1}, {'Xlim', 1}, {'Zlim', 1}, {'NoPause'}}, varargin);
            
            equalz = opts(1);
            if (opts(2))
                titl = args{2};
            else
                titl = -1;
            end
            
            if (opts(3))
                ylab = args{3};
            else
                ylab = -1;
            end
            
            if (opts(4))
                xlimin = true;
                xli = args{4};
            else
                xlimin = false;
                xli = [];
            end
            
            if (opts(5))
                zlimin = true;
                zli = args{5};
            else
                zlimin = false;
                zli = [];
            end
            
            noPause = opts(6);
                        
            ti = 0:dt:runTime;
            Nti = length(ti);
            
            Nmvs = length(motionsVids);
            
            figure;
            
            limx = zeros(Nmvs, 1);
            limz = zeros(Nmvs, 1);
            maxz = 0;
            for n = 1:Nmvs
                %plotaxes(n) = subplot(Nmvs, 1, n);
                
                if (xlimin)
                    limx(n) = xli;
                    motionsVids(n).Xlim = xli;
                else
                    limx(n) = motionsVids(n).Xlim;
                end
                
                if (zlimin)
                    limz(n) = zli;
                else
                    limz(n) = motionsVids(n).Zlim;
                end
                
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
            
            thisFig = gcf;
            for m = 1:Nti
                for n = 1:Nmvs
                    if (gcf ~= thisFig)
                        return;
                    end
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
                        plot(xwav,wav, 'LineWidth', 1.2);
                    end
                    hold on;
                    
                    if (iscell(ylab))
                        ylabel(ylab{n});
                    end
                    
                    if (iscell(bp))
                        for o = 1:length(bp)
                            scatter(bp{o}(:,1), bp{o}(:,3), '.');
                            plot(bp{o}(:,1), bp{o}(:,3), 'LineWidth', 1.2 );
                        end
                    else
                        scatter(bp(:,1), bp(:,3), '.');
                        plot(bp(:,1), bp(:,3), 'LineWidth', 1.2 );
                    end
                    if (m == 1)
                        axis equal
                    end
                end
                if (~noPause)
                    pause(0.5*dt);
                end
                mov(m) = getframe(gcf);
                
            end
        end
    end
    
    methods (Access = private)
        function [] = computeMotWave(mv)            
            mv.omega = 2*pi/mv.t;                      
                        
            if (mv.nBod == 1)
                Np = size(mv.bodyPoints, 1);

                mv.pointMot = zeros(Np, 3);

                for m = 1:Np
                    for n = 1:mv.dof
                        mv.pointMot(m,:) = mv.pointMot(m,:) + mv.motions(n)*mv.motionFuncs(n).Evaluate(mv.bodyPoints(m,:));
                    end
                end

                bp = mv.bodyPoints + abs(mv.pointMot);
            else
                mv.pointMot = cell(mv.nBod, 1);
                bp = cell(mv.nBod, 1);
                for l = 1:mv.nBod
                    Np = size(mv.bodyPoints{l}, 1);

                    mv.pointMot{l} = zeros(Np, 3);

                    for m = 1:Np
                        for n = 1:mv.dof
                            mv.pointMot{l}(m,:) = mv.pointMot{l}(m,:) + mv.motions{l}(n)*mv.motionFuncs{l}(n).Evaluate(mv.bodyPoints{l}(m,:));
                        end
                    end

                    bp{l} = mv.bodyPoints{l} + abs(mv.pointMot{l});
                end
            end
            
            k = IWaves.SolveForK(mv.omega, mv.h, IWaves.G);
            
            if (isempty(mv.xlim))
                lam = 2*pi./k;
                if (mv.nBod == 1)
                    limx = max(abs(bp(:,1)));
                else
                    limx = 0;
                    for n = 1:mv.nBod
                        limxn = max(abs(bp{n}(:,1)));
                        if (limxn > limx)
                            limx = limxn;
                        end
                    end
                end
                if (limx < 0.25*lam)
                    limx = 0.25*lam;
                end
                mv.xlim = limx + 0.5*limx;
            end
            
            if (mv.nBod == 1)
                limz = max(abs(bp(:,3)));
            else
                limz = 0;
                for n = 1:mv.nBod
                    limzn = max(abs(bp{n}(:,3)));
                    if (limzn > limz)
                        limz = limzn;
                    end
                end
            end
%             if (limz < 1)
%                 limz = 1;
%             end
            
            mv.zlim = limz + 0.5*limz;
            
            if (mv.showWave)
                if (isempty(mv.h))
                    error('Depth value (H) must be set to show the wave');
                end
               
                mv.xWave = -mv.xlim:0.1:mv.xlim;
                mv.wave = mv.waveA*exp(-1i*k*mv.xWave);
            end
            
            mv.isComp = true;
        end
    end
end