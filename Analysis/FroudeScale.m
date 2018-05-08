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
classdef FroudeScale < handle    
    methods (Static)
        function [vals] = Compute(type, scale, val)
            scl = [];
            
            if strcmpi(type, 'linear')
                scl = scale;
            elseif strcmpi(type, 'area')
                scl = scale.^2;
            elseif strcmpi(type, 'volume')
                scl = scale.^3;
            elseif strcmpi(type, 'angle')
                scl = 1;
            elseif strcmpi(type, 'time')
                scl = sqrt(scale);
            elseif strcmpi(type, 'frequency')
                scl = 1./sqrt(scale);
            elseif strcmpi(type, 'vel')
                scl = sqrt(scale);
            elseif strcmpi(type, 'angvel')
                scl = 1./sqrt(scale);
            elseif strcmpi(type, 'accel')
                scl = 1;
            elseif strcmpi(type, 'angaccel')
                scl = 1./sqrt(scale);
            elseif strcmpi(type, 'mass')
                scl = scale.^3;
            elseif strcmpi(type, 'inertia')
                scl = scale.^5;
            elseif strcmpi(type, 'stifflin')
                scl = scale.^2;
            elseif strcmpi(type, 'stiffrot')
                scl = scale.^4;
            elseif strcmpi(type, 'damplin')
                scl = scale.^2.5;
            elseif strcmpi(type, 'damprot')
                scl = scale.^4.5;
            elseif strcmpi(type, 'force')
                scl = scale.^3;
            elseif strcmpi(type, 'torque')
                scl = scale.^4;
            elseif strcmpi(type, 'power')
                scl = scale.^3.5;
            elseif strcmpi(type, 'energy')
                scl = scale.^4;
            else
                error('scale type not found');
            end
            
            vals = scl*val;
        end
    end
end