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
classdef ConstraintMatComp
    
    methods (Static)
        function [P] = HingedBodies(bods, hins, varargin)
            % bods = N x 3 matrix of {x, y, z} body coordinates, where N is
            % the number of bodies
            % hins = (N-1) x 3 matrix of {x, y, z} hinge coordinates. The
            % y-coordinate is not really necessary as it's a hinge about an
            % axis parallel to the y-axis
            
            [opts, args] = checkOptions({{'Origin', 1}}, varargin);
            
            if (opts(1))
                org = args{1};
            else
                org = bods(1,:);
            end
            
            Nbod = size(bods, 1);
            Nhin = size(hins, 1);
            
            if (Nhin ~= (Nbod - 1))
                error('The number of hinges must be one less than the number of bodies');
            end
            
            if (Nbod > 1)
                if ((size(bods, 2) ~= 3) || (size(hins, 2) ~= 3))
                    error('The body coordinates and hinge coordinates must have x,y,z locations');
                end
            end
            
            PT = zeros(6*Nbod, 6 + Nhin);
            
            sR = zeros(Nbod, 3);
            sL = zeros(Nbod, 3);
            
            for n = 1:Nbod
                if (n < Nbod)
                    sR(n,:) = hins(n,:) - bods(n,:);
                end
                
                if (n > 1)
                    sL(n,:) = hins(n-1,:) - bods(n,:);
                end
            end
            
            sO = bods(1,:) - org;
              
            for n = 1:Nbod
                PTn = zeros(6, 6 + Nhin);
                for m = 1:6
                    PTn(m, m) = 1;
                end
                
                svect = -sO;
                
                for m = 2:n
                    svect = svect - sR(n-1,:) + sL(n,:);
                end
                
                Sx = ConstraintMatComp.skewMat(svect);
                PTn(1:3,4:6) = Sx;
                
                PTn(5,7:(5+n)) = ones(1,n-1);
                
                if (n > 1)
                    svect = [-sL(n,3), 0, sL(n,1)];
                    for m = (n-1):-1:2
                        svect = svect + [-sL(m,3), 0, sL(m,1)] - [-sR(m,3), 0, sR(m,1)];
                    end
                    PTn(1:3,n+5) = svect';
                end
                
                istart = (n - 1)*6 + 1;
                PT(istart:(istart + 5), :) = PTn;
            end
            
            P = PT.';
        end
    end
    
    methods (Static, Access = private)
        
        function [M] = skewMat(v)
            M = zeros(3, 3);
            
            M(1,2) = -v(3);
            M(1,3) = v(2);
            M(2,1) = v(3);
            M(2,3) = -v(1);
            M(3,1) = -v(2);
            M(3,2) = v(1);
        end
    end
end