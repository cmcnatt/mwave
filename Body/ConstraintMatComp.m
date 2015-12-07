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
        function [P] = Hinge(bods, hins, varargin)
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
            
            if ((size(bods, 2) ~= 3) || (size(hins, 2) ~= 3))
                error('The body coordinates and hinge coordinates must have x,y,z locations');
            end
            
            PT = zeros(6*Nbod, 6 + Nhin);
            
            dhbR = zeros(3, Nbod);
            dhbL = zeros(3, Nbod);
            
            for n = 1:Nbod
                if (n < Nbod)
                    dhbR(n,:) = hins(n,:) - bods(n,:);
                end
                
                if (n > 1)
                    dhbL(n,:) = hins(n-1,:) - bods(n,:);
                end
            end
            
            dorgh1 = hin(1,:) - org;
                        
            for n = 1:Nbod
                PTn = zeros(6, 6 + Nhin);
                for m = 1:6
                    PTn(m, m) = 1;
                end
                
                PTn(1,5) = dorgh1(3);
                PTn(1,6) = dorgh1(2);
                PTn(2,4) = dorgh1(3);
                PTn(2,6) = dorgh1(1);
                PTn(3,4) = dorgh1(2);
                PTn(3,5) = dorgh1(1);
                
                for m = 1:n
                    PTn(1,5) = PTn(1,5) + dhbL(m,3) - dhbR(m,3);
                    PTn(1,6) = PTn(1,6) + dhbL(m,2) - dhbR(m,2);
                    PTn(2,4) = PTn(2,4) + dhbL(m,3) - dhbR(m,3);
                    PTn(2,6) = PTn(2,6) + dhbL(m,1) - dhbR(m,1);
                    PTn(3,4) = PTn(3,4) + dhbL(m,2) - dhbR(m,2);
                    PTn(3,5) = PTn(3,5) + dhbL(m,1) - dhbR(m,1);
                end
                
                PTn(1,6) = -PTn(1,6);
                PTn(2,4) = -PTn(2,4);
                PTn(3,5) = -PTn(3,5);
                
                for m = (6+Nhin):n
                    PTn(5, m) = 1;
                end
                                
                istart = (n - 1)*6 + 1;
                PT(istart:(istart + 6)) = PTn;
            end
        end
    end
end