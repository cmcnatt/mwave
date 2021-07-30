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
classdef StlVideo < handle
    
    properties (Access = private)
        bods;
        mots;
        wave;
        time;
        ramp;
    end
    
    properties (Dependent)
        Bodies;
        Motions;
        Wave;
        Time;
        Ramp;
    end
    
    methods
        function [vid] = StlVideo()
        end
        
        function [] = set.Bodies(stl, val)
            stl.bods = val;
        end
        function [val] = get.Bodies(stl)
            val = stl.bods;
        end
        
        function [] = set.Motions(stl, val)
            stl.mots = val;
        end
        function [val] = get.Motions(stl)
            val = stl.mots;
        end
        
        function [] = set.Wave(stl, val)
            stl.wave = val;
        end
        function [val] = get.Wave(stl)
            val = stl.wave;
        end
        
        function [] = set.Time(stl, val)
            stl.time = val;
        end
        function [val] = get.Time(stl)
            val = stl.time;
        end
        
        function [] = set.Ramp(stl, val)
            stl.ramp = val;
        end
        function [val] = get.Ramp(stl)
            val = stl.ramp;
        end
        
        function [mov] = movie(vid, varargin)
            % overloaded movie function
            
            [opts, args] = checkOptions({{'pause'}, {'view', 1}, {'write', 1}, {'vids', 1}}, varargin);
            noPause = ~opts(1);
            ang = [0, 0];
            if opts(2)
                ang = args{2};
            end
            write = false;
            if opts(3)
                write = true;
                vidFile = args{3};
            end
            otherVids = {};
            if opts(4)
                otherVids = args{4};
            end
            
            Nvid = 1;
            if ~isempty(otherVids)
                Nvid = 1 + length(otherVids);
            end
                                    
            Nti = length(vid.time);
            limx = [min(vid.wave.X) max(vid.wave.X)];
            limy = [min(vid.wave.Y) max(vid.wave.Y)];
            
            for n = 1:length(otherVids)
                vidn = otherVids{n};
                if length(vid.Time) ~= Nti
                    error('All videos must have same time length');
                end
                limxn = [min(vidn.Wave.X) max(vidn.Wave.X)];
                limyn = [min(vidn.Wave.Y) max(vidn.Wave.Y)];
                if ~all(limxn == limx) && ~all(limyn == limy)
                    error('All videos must have same size wave');
                end
            end
            
            vids = {vid, otherVids{:}};
            Nvid = length(vids);
            
            maxz = 0;
            minz = 0;
            
            for m = 1:length(vids)
                bodsm = vids{m}.Bodies;
                if ~iscell(bodsm)
                    bodsm = {bodsm};
                end

                for n = 1:length(bodsm)
                    bodsm{n}.SetSetting(zeros(1, bodsm{n}.DOF));
                    [~, cents, ~] = bodsm{n}.PointsAtSetting();
                    if iscell(cents)
                        z = cents{1}(:,3);
                        for m = 2:length(cents)
                            zm = cents{m}(:,3);
                            z(length(z)+1:length(z)+length(zm)) = zm; 
                        end
                    else
                        z = cents(:,3);
                    end
                    if max(z) > maxz
                        maxz = max(z);
                    end
                    if min(z) < minz
                        minz = min(z);
                    end
                end
            end
            
            limz = [2*minz 2*maxz];
            
            [maxz, minz] = vid.wave.ComputeMaxMin(vid.time);
            limw = [minz maxz];
            
            mov = [];
            nFrames = Nti;
            if write
                v = VideoWriter(vidFile);
                v.FrameRate = 1/(vid.time(2) - vid.time(1));
                open(v);
            else
%                 mov(1:nFrames) = struct('cdata',[], 'colormap',[]);
                for i = 1:nFrames; mov(i).cdata = [];mov(i).colormap = []; end; % this line was switched on 30-7-21 - it seems as if a Matlab update has prevented the format used previously.
            end
            
            figure;
            thisFig = gcf;
            
            for ti = 1:Nti
                if (gcf ~= thisFig)
                    close(v);
                    return;
                end
                for m = 1:Nvid
                    vidm = vids{m};
                    subplot(Nvid, 1, m);
                    cla;
                    axis equal;
                    set(gca, 'xlim', limx, 'ylim', limy, 'zlim', limz, 'clim', limw);
                    if ti == 1
                        set(gca, 'view', ang);
                    end

                    if ~isempty(vidm.Ramp)
                        vidm.Wave.Ramp = vidm.Ramp(ti);
                    end
                    vidm.Wave.Time = vidm.Time(ti);
                    plot(vidm.Wave);
                    
                    bodsm = vidm.Bodies;
                    if ~iscell(bodsm)
                        bodsm = {bodsm};
                    end
                    prevDOF = 0; % Initialise DOF counter
                    for n = 1:length(bodsm)
                        bodsm{n}.SetSetting(vidm.Motions(ti,1+(n-1)*bodsm{n}.DOF:n*bodsm{n}.DOF));
                        bodsm{n}.SetSetting(vidm.Motions(ti,prevDOF+1:prevDOF+bodsm{n}.DOF));
                        prevDOF = prevDOF+bodsm{n}.DOF; % This should allow for bodies with any amounts of DOFs.
                        plot(bodsm{n}, 'Color', MColor.Black);
                    end
                end
                if (~noPause)
                    pause(0.5*dt);
                end
                if write
                    writeVideo(v, getframe(gcf));
                else
                    mov(ti) = getframe(gcf);
                end
            end
            
            if write
                close(v);
            end
        end
    end
end