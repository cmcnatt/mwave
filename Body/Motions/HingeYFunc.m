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
classdef HingeYFunc < IMotionFunc
    % Defines motions for modes of motion
    properties (Access = private)
        hingePos;
    end
    
    properties (Dependent)
        MotionIn;
        HingePos;
    end
    
    methods
        function [hinge] = HingeYFunc(varargin)
            hinge.initCg(varargin{:});
            hinge.hingePos = [0 0];
            hinge.isym = 1;
            hinge.isGen = true;
        end
        
        function [f] = Evaluate(hinge, pos)
            hinPos = hinge.hingePos;
            dx = pos(1) - hinPos(1);
            fx = abs(pos(3) - hinPos(2))*sign(dx);
            fy = 0;
            fz = abs(dx);
            
            f = [fx fy fz];
        end
        
        function [div] = Divergence(hinge, pos)
            div = 0;
        end
        
        function [gf] = GravityForce(hinge, pos)
            cg_ = hinge.cg;
            pos = pos - cg_;
            pos = [abs(pos(1)) pos(2:3)];
            znorm = [0 0 -1];
            gf = cross(pos, znorm);
            gf = gf(2);
        end
        
        function [in] = get.MotionIn(hinge)
            in = [1 0 1];
        end
        
        function [hp] = get.HingePos(hinge)
            hp = hinge.hingePos;
        end
        function [hinge] = set.HingePos(hinge, hp)
            [row, col] = size(hp);
            if ((row ~= 1) || (col ~= 2))
                error('HingePos must be a 1x2 vector of the {x,z} location of the hinge');
            end
            
            hinge.hingePos = hp;
        end
    end
end