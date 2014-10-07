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
classdef ZeroWaveField < WaveField
    
    methods
        function [wf] = ZeroWaveField(varargin)
            if (length(varargin) == 1)
                wfin = varargin{1};
                rho = wfin.Rho;
                h = wfin.H;
                t = wfin.T;
                isarray = wfin.IsArray;
                
                if (isarray)
                    [X, Y] = wfin.FieldPoints;
                    varargin = {X, Y};
                else
                    pts = wfin.FieldPoints;
                    varargin = {pts};
                end
                
                if(~wfin.HasVelocity)
                    noVel = true;
                else
                    noVel = false;
                end
            else
                rho = varargin{1};
                h = varargin{2};
                t = varargin{3};
                isarray = varargin{4};
                varargin = varargin(5:end);
                
                noVel = checkOptions({'NoVel'}, varargin);
            end
            
            nT = length(t);
            [nY, nX] = size(varargin{1});
            if (isarray)
                p = zeros(nT, nY, nX);
            else
                p = zeros(nT, nY);
            end
            
            if (noVel)
                vel = [];
            else
                if (isarray)
                    vel = zeros(nT, 3, nY, nX);
                else
                    vel = zeros(nT, 3, nY);
                end
            end
            
            wf = wf@WaveField(rho, IWaves.G, h, t, p, vel, isarray, varargin{:});
        end
    end
    
end