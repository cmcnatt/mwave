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
classdef ViscDampCoef < handle
    % Defines the damping coefficient for a Morison type body-motion
    % viscous damping
    %
    % F = 1/2*rho*cd*A*|u|u
    %
    % rho - fluid density
    % cd - non-dimensional damping coefficient
    % A - characteristic "area", for rotations, it has dimensions of L^5
    % u - body velocity
    
    properties (Access = private)
        cd;
        rho;
        a;
        isRot;
    end
    
    properties (Dependent)
        Cd;                 % The dimensionless damping coefficient
        Rho;                % The fluid density
        A;                  % The area (trans), the  third moment of the area (dim: len^5) (rotational)
        IsRotational;       % Indicates whether the Damping is rotational
    end
    
    methods
        
        function [coef] = ViscDampCoef()
            
            coef.rho = [];
            coef.cd = 0;
            coef.a = 0;
            coef.isRot = [];
        end
        
        function [val] = get.Cd(coef)
            % The dimensionless damping coefficient
            val = coef.cd;
        end
        function [] = set.Cd(coef, val)
            % The dimensionless damping coefficient
            if double(val) < 0 || ~isscalar(double(val))
                error('Cd must be a positive scalar');
            end
            coef.cd = val;
        end
        
        function [val] = get.Rho(coef)
            % The fluid density
            val = coef.rho;
        end
        function [] = set.Rho(coef, val)
            % The fluid density
            if val < 0 || ~isscalar(val)
                error('Rho must be a positive scalar');
            end
            coef.rho = val;
        end
        
        function [val] = get.A(coef)
            % The area (trans), the  third moment of the area (dim: len^5) (rotational)
            val = coef.a;
        end
        function [] = set.A(coef, val)
            % The area (trans), the  third moment of the area (dim: len^5) (rotational)
            if double(val) < 0 || ~isscalar(double(val))
                error('A must be a positive scalar');
            end
            coef.a = val;
        end
        
        function [val] = get.IsRotational(coef)
            % Indicates whether the Damping is rotational
            val = coef.isRot;
        end
        function [] = set.IsRotational(coef, val)
            % Indicates whether the Damping is rotational
            if ~isBool(val) || ~isscalar(val)
                error('IsRotational must be a scalar boolean');
            end
            coef.isRot = val;
        end
        
        function [val] = GetDimCd(coef, varargin)
            % The dimensional damping coefficient:
            %   CD = 1/2*rho*cd*A
            % The Lorentz linearisation (Ref: Folley & Whittaker (2010),
            % Spectral modelling of wave energy converters. Eq. 13)
            %   CD = (8/3*pi)*1/2*rho*cd*A
            
            opts = checkOptions({{'linear'}}, varargin);
            lin = opts(1);
            
            cd_ = IRandomVariable.TrySample(coef.cd, varargin{:});
            a_ = IRandomVariable.TrySample(coef.a, varargin{:});
            
            val = 1/2*coef.rho*cd_.*a_;
            if lin
                val = 8/(3*pi)*val;
            end
        end
    end
    
    methods (Static)
        function [A] = ComputeAreaFromBox(x, y, z, centRot, varargin)
            % assumes the area used in the damping coefficient calculation
            % can be computed from the bounding box around the geometry
            %
            % The translational areas are just the projected box areas in
            % the plane orthogonal to the direction
            %
            % For the rotational area:
            %  - we need dimentions of length^5
            %  - the velocity at a distance r is: u = r*theta_dot, where
            % theta_dot is the rotational velocity, and r is the distance
            % from the center of rotation
            % - the viscous damping is proportional to |u|u, so that means
            % for rotation it is proportional to |r|r
            % - it acts over an area L*dr, where L is the orthogonal length,
            % and dr is an infintesimal area
            % - it creats a torque: T = F*r
            % - Finally, we have the integral:
            %      
            %   A = L?|r|^3dr
            %     = 1/4*L*r^4*sgn(r)
            
            [opts, args] = checkOptions({{'makeCoefs'}, {'rho', 1}}, varargin);
            makeCoefs = opts(1);
            rho = [];
            if opts(2)
                rho = args{2};
            end
            
            dx = x(2) - x(1);
            dy = y(2) - y(1);
            dz = z(2) - z(1);
            
            rx = [-sqrt((y(1)-centRot(2))^2 + (z(1)-centRot(3))^2), ...
                sqrt((y(2)-centRot(2))^2 + (z(2)-centRot(3))^2)];
            
            ry = [-sqrt((x(1)-centRot(1))^2 + (z(1)-centRot(3))^2), ...
                sqrt((x(2)-centRot(1))^2 + (z(2)-centRot(3))^2)];
            
            rz = [-sqrt((x(1)-centRot(1))^2 + (y(1)-centRot(2))^2), ...
                sqrt((x(2)-centRot(1))^2 + (y(2)-centRot(2))^2)];
            
            A = zeros(1, 6);
            
            A(1) = dy*dz;
            A(2) = dx*dz;
            A(3) = dx*dy;
            A(4) = dx/4*(rx(2)^4+rx(1)^4);
            A(5) = dy/4*(ry(2)^4+ry(1)^4);
            A(6) = dz/4*(rz(2)^4+rz(1)^4);
            
            if makeCoefs
                coefs(1, 6) = ViscDampCoef;
                for n = 1:6
                    coefs(n).A = A(n);
                    if n <= 3
                        coefs(n).IsRotational = false;
                    else
                        coefs(n).IsRotational = true;
                    end
                    if ~isempty(rho)
                        coefs(n).Rho = rho;
                    end
                end
                A = coefs;
            end
        end
    end
end