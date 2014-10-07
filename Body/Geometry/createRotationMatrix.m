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
function [rotMat] = createRotationMatrix(ax, angle)

cAng = cos(angle);
sAng = sin(angle);

lenAxis = sqrt(ax(1)^2 + ax(2)^2 + ax(3)^2);

ux = ax(1)/lenAxis;
uy = ax(2)/lenAxis;
uz = ax(3)/lenAxis;

rotMat = zeros(3,3);

rotMat(1,1) = cAng + ux^2*(1-cAng);
rotMat(1,2) = ux*uy*(1-cAng) - uz*sAng;
rotMat(1,3) = ux*uz*(1-cAng) + uy*sAng;

rotMat(2,1) = uy*ux*(1-cAng) + uz*sAng;
rotMat(2,2) = cAng + uy^2*(1-cAng);
rotMat(2,3) = uy*uz*(1-cAng) - ux*sAng;

rotMat(3,1) = uz*ux*(1-cAng) - uy*sAng;
rotMat(3,2) = uz*uy*(1-cAng) + ux*sAng;
rotMat(3,3) = cAng + uz^2*(1-cAng);

end
