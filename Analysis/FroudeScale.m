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
            elseif strcmpi(type, 'addedmasslinlin')
                scl = scale.^3;
            elseif strcmpi(type, 'addedmasslinrot')
                scl = scale.^4;
            elseif strcmpi(type, 'addedmassrotrot')
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
            elseif strcmpi(type, 'forcerao')
                scl = scale.^2;
            elseif strcmpi(type, 'torque')
                scl = scale.^4;
            elseif strcmpi(type, 'torquerao')
                scl = scale.^3;
            elseif strcmpi(type, 'power')
                scl = scale.^3.5;
            elseif strcmpi(type, 'powerrao')
                scl = scale.^1.5;
            elseif strcmpi(type, 'energy')
                scl = scale.^4;
            else
                error('scale type not found');
            end
            
            vals = scl*val;
        end
        
        function [vals] = Linear(scale, val)
            vals = FroudeScale.Compute('linear', scale, val);
        end
        
        function [vals] = Area(scale, val)
            vals = FroudeScale.Compute('area', scale, val);
        end
        
        function [vals] = Volume(scale, val)
            vals = FroudeScale.Compute('volume', scale, val);
        end
        
        function [vals] = Angle(scale, val)
            vals = FroudeScale.Compute('angle', scale, val);
        end
        
        function [vals] = Time(scale, val)
            vals = FroudeScale.Compute('time', scale, val);
        end
        
        function [vals] = Velocity(scale, val)
            vals = FroudeScale.Compute('vel', scale, val);
        end
        
        function [vals] = AngularVelocity(scale, val)
            vals = FroudeScale.Compute('angvel', scale, val);
        end
        
        function [vals] = Acceleration(scale, val)
            vals = FroudeScale.Compute('accel', scale, val);
        end
        
        function [vals] = AngularAcceleration(scale, val)
            vals = FroudeScale.Compute('angaccel', scale, val);
        end
        
        function [vals] = Mass(scale, val)
            vals = FroudeScale.Compute('mass', scale, val);
        end
        
        function [vals] = Inertia(scale, val)
            vals = FroudeScale.Compute('inertia', scale, val);
        end
        
        function [vals] = StiffnessLinear(scale, val)
            vals = FroudeScale.Compute('stifflin', scale, val);
        end
        
        function [vals] = StiffnessRotational(scale, val)
            vals = FroudeScale.Compute('stiffrot', scale, val);
        end
        
        function [vals] = DampingLinear(scale, val)
            vals = FroudeScale.Compute('damplin', scale, val);
        end
        
        function [vals] = DampingRotational(scale, val)
            vals = FroudeScale.Compute('damprot', scale, val);
        end
        
        function [vals] = Force(scale, val)
            vals = FroudeScale.Compute('force', scale, val);
        end
        
        function [vals] = ForceRAO(scale, val)
            vals = FroudeScale.Compute('forcerao', scale, val);
        end
        
        function [vals] = Torque(scale, val)
            vals = FroudeScale.Compute('torque', scale, val);
        end
        
        function [vals] = TorqueRAO(scale, val)
            vals = FroudeScale.Compute('torquerao', scale, val);
        end
        
        function [vals] = Power(scale, val)
            vals = FroudeScale.Compute('power', scale, val);
        end
        
        function [vals] = PowerRAO(scale, val)
            vals = FroudeScale.Compute('powerrao', scale, val);
        end
        
        function [vals] = Energy(scale, val)
            vals = FroudeScale.Compute('energy', scale, val);
        end
    end
end