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
    end
    
    properties (Dependent)
        Bodies;
        Motions;
        Wave;
        Time;
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
        
        function [mov] = movie(vid, varargin)
            % overloaded movie function
            
            [opts, args] = checkOptions({{'pause'}, {'view', 1}}, varargin);
            noPause = ~opts(1);
            ang = [0, 0];
            if opts(2)
                ang = args{2};
            end
            
            Nti = length(vid.time);
            
            bods_ = vid.bods;
            if ~iscell(bods_)
                bods_ = {bods_};
            end
            
            figure;
                        
            limx = [min(vid.wave.X) max(vid.wave.X)];
            limy = [min(vid.wave.Y) max(vid.wave.Y)];
            
            maxz = 0;
            minz = 0;
            for n = 1:length(bods_)
                bods_{n}.SetSetting(zeros(1, bods_{n}.DOF));
                [~, cents, ~] = bods_{n}.PointsAtSetting();
                z = cents(:,3);
                if max(z) > maxz
                    maxz = max(z);
                end
                if min(z) < minz
                    minz = min(z);
                end
            end
            
            limz = [3*minz 3*maxz];
            
            nFrames = Nti;
            mov(1:nFrames) = struct('cdata',[], 'colormap',[]);
            
            thisFig = gcf;
            
            for m = 1:Nti
                if (gcf ~= thisFig)
                    return;
                end
                cla;
                axis equal;
                set(gca, 'xlim', limx, 'ylim', limy, 'zlim', limz, 'clim', limz);
                if m == 1
                    set(gca, 'view', ang);
                end
                
                vid.wave.Time = vid.time(m);
                plot(vid.wave);
                for n = 1:length(bods_)
                    bods_{n}.SetSetting(vid.mots(m,:));
                    plot(bods_{n});
                end
                if (~noPause)
                    pause(0.5*dt);
                end
                mov(m) = getframe(gcf);
                
            end
        end
    end
end