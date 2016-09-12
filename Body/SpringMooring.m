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
classdef SpringMooring < handle
    
    properties (Access = protected)
        loc;
        k;
        body;
        angsHor;
        angsVer;
    end
    
    properties (Dependent)
        Location;
        LineCount;
        SpringConst;
        BodyAttached;
        AngsHor;
        AngsVer;
        MooringMatrix;
    end
    
    methods 
        
        function [moor] = SpringMooring(varargin)
            if ~isempty(varargin)
                moor.loc = varargin{1};
                moor.k = varargin{2};
                moor.angsHor = varargin{3};
                moor.angsVer = varargin{4};
                if length(varargin) > 4
                    moor.body = varargin{5};                    
                    moor.SetMooringMatrix;
                end
            end
        end
        
        function [] = set.AngsHor(moor, val)
            moor.angsHor = val;
        end
        function [val] = get.AngsHor(moor)
            val = moor.angsHor;
        end
        
        function [] = set.AngsVer(moor, val)
            moor.angsVer = val;
        end
        function [val] = get.AngsVer(moor)
            val = moor.angsVer;
        end
        
        function [] = set.BodyAttached(moor, val)
            moor.body = val;
        end
        function [val] = get.BodyAttached(moor)
            val = moor.body;
        end
        
        function [] = set.Location(moor, val)
            moor.loc = val;
        end
        function [val] = get.Location(moor)
            val = moor.loc;
        end
        
        function [val] = get.LineCount(moor)
            val = length(moor.angsHor);
        end
        
        function [] = set.SpringConst(moor, val)
            moor.k = val;
        end
        function [val] = get.SpringConst(moor)
            val = moor.k;
        end
        
        function [K] = SetMooringMatrix(moor)
            N = length(moor.angsHor);
            if 0 == N
                error('Not all properties are set in the SpringMooring to compute to mooring matrix');
            end
            if length(moor.angsVer) ~= N
                error('The number of horizonatal angles must equal the number of vertical angles');
            end
            if isempty(moor.body)
                error('The AttachedBody is not set');
            end
            if isempty(moor.loc)
                error('The Location is not set');
            end
            if isempty(moor.k)
                error('The spring constant is not set');
            end
            
            Cg = moor.body.Cg;
            kpos = moor.loc;
            
            if 1 == length(moor.k)
                ks = moor.k*[1, 1, 1];      % N/m
            elseif N == length(moor.k)
                ks = moor.k;
            else
                error('The number of spring constants must be 1 or equal to the number of mooring lines');
            end
            
            anghs = moor.angsHor;
            angvs = moor.angsVer;
            
            K = computeMooringK(Cg, kpos, ks, anghs, angvs);
            
            moor.body.K = K;
        end
        
        function [] = PlotMooring(moor, ax, depth, is2d, varargin)
            N = length(moor.angsHor);
            pt1 = moor.loc;
            if isinf(depth)
                depth = -pt1(3)+10;
            end
            
            z = -depth;
            for n = 1:N
                linLen = (depth+pt1(3))/sin(moor.angsVer(n));
                x = pt1(1) + linLen*cos(moor.angsHor(n));
                y = pt1(2) + linLen*sin(moor.angsHor(n));
                hold on;
                if is2d
                    plot(ax, [pt1(1) x], [pt1(3) z], varargin{:});
                else
                    plot3(ax, [pt1(1) x], [pt1(2) y], [pt1(3) z], varargin{:});
                end
            end
        end
    end
    
    
end