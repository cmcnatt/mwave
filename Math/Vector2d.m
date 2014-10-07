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
classdef Vector2d
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        x;
        y;
    end
    
    properties (Dependent)
        X
        Y
        R
        Theta
    end
    
    methods
        function [vect] = Vector2d(varargin)
            if (~isempty(varargin))
                if (length(varargin) ~= 2)
                    error('there must be two inputs');
                end
                if (~isscalar(varargin{1}))
                    error('input must be a scalar value');
                else
                    vect.x = varargin{1};
                end
                if (~isscalar(varargin{2}))
                    error('input must be a scalar value');
                else
                    vect.y = varargin{2};
                end
            else
                vect.x = 1;
                vect.y = 0;
            end
        end
        
        function [x_] = get.X(vect)
            x_ = vect.x;         
        end
        function [vect] = set.X(vect, x_)
            if (~isscalar(x_))
                error('input must be a scalar value');
            end
            vect.x = x_;
        end
        
        function [y_] = get.Y(vect)
            y_ = vect.y;         
        end
        function [vect] = set.Y(vect, y_)
            if (~isscalar(y_))
                error('input must be a scalar value');
            end
            vect.y = y_;
        end
        
        function [r] = get.R(vect)
            r = sqrt(vect.x^2 + vect.y^2);      
        end
        function [vect] = set.R(vect, r)
            if (~isscalar(r))
                error('input must be a scalar value');
            end
            thet = vect.Theta;
            vect.x = r*cos(thet);
            vect.y = r*sin(thet);
        end
        
        function [theta] = get.Theta(vect)
            theta = atan2(vect.y, vect.x);    
        end
        function [vect] = set.Theta(vect, theta)
            if (~isscalar(theta))
                error('input must be a scalar value');
            end
            r = vect.R;
            vect.x = r*cos(theta);
            vect.y = r*sin(theta);
        end
        
        function [areEq] = eq(va, vb)
            if (~isa(va, 'Vector2d') || ~isa(vb, 'Vector2d'))
                error('Each argument must be a VEctor2d');
            end
            
            if (va.x == vb.x && va.y == vb.y)
                areEq = true;
            else
                areEq = false;
            end
        end
        
        function [areNE] = ne(va, vb)
            if (va == vb)
                areNE = false;
            else
                areEq = true;
            end
        end
        
        function [vout] = plus(va, vb)
            if (~isa(va, 'Vector2d') || ~isa(vb, 'Vector2d'))
                error('Each argument must be a VEctor2d');
            end
            
            vout = Vector2d(va.x + vb.x, va.y + vb.y);
        end
        
        function [vout] = uminus(vin)
            if (~isa(vin, 'Vector2d'))
                error('Argument must be a VEctor2d');
            end
            
            vout = Vector2d(-vin.x, -vin.y);
        end
        
        function [vout] = minus(va, vb)
            vout = va + (-vb);
        end
        
        function [out] = mtimes(va, vb)
            scal = false;
            if (~isa(va, 'Vector2d'))
                if (isscalar(va))
                    error('Argument must be a Vector2d or a scalar');
                else
                    scal = true;
                    vscal = va;
                    vvect = vb;
                end
            end
            
            if (~isa(vb, 'Vector2d'))
                if (isscalar(vb))
                    error('Argument must be a Vector2d or a scalar');
                else
                    scal = true;
                    vscal = vb;
                    vvect = va;
                end
            end
            
            if (scal)
                out = Vector2d(vscal*vvect.x, vscal*vvect.y);
            else
                out = va.x*vb.x + va.y*vb.y;
            end
        end
        
        function [out] = times(va, vb)
            out = va*vb;
        end
    end
    
end

