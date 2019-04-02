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
classdef ForceComp
    
    methods (Static)
        function [Fout] = ForceP2O(Fin, point, origin)
            if nargin < 3
                origin = [0; 0; 0];
            end
            if isrow(origin)
                origin = origin';
            end
            if isrow(Fin)
                Fin = Fin.';
            end
            if length(Fin) == 3
                Fin = [Fin; zeros(3,1)];
            end
            if isrow(point)
                point = point.';
            end
            f = Fin(1:3);
            r = point - origin;
            
            rskew = skewMat(r);
            
            Fout = zeros(6, 1);
            Fout(1:3) = f;
            Fout(4:6) = Fin(4:6) + cross(r, f);
        end
    end
end