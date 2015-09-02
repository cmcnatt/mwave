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
classdef YawFunc < IMotionFunc
    % Defines motions for modes of motion
    
    properties (Dependent)
        MotionIn;
    end
    
    methods
        function [yaw] = YawFunc(varargin)
            yaw.initCg(varargin{:});
            yaw.isym = 6;
            yaw.isGen = false;
        end
        
        function [f] = Evaluate(yaw, pos)
            if (yaw.checkEvalPoint(pos))
                cg_ = yaw.cg;
                fx = -(pos(2) - cg_(2));
                fy = pos(1) - cg_(1);
                fz = 0;

                f = [fx fy fz];
            else
                f = [0 0 0];
            end            
        end
        
        function [div] = Divergence(yaw, pos)
            div = 0;
        end
        
        function [gf] = GravityForce(yaw, pos)
            gf = 0;
        end
        
        function [in] = get.MotionIn(yaw)
            in = [1 1 0];
        end
    end
end