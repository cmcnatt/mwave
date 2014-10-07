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
function [r0, theta, z, eta] = reshapeCirWFPoints(points, eta0)
% Reshapes the listed circular point locations and wave elevations into a 
% rectangular matrix.  Needed for wave decomposition (see waveDecomp)

r = sqrt(points(:,1).^2 + points(:,2).^2);

r0 = r(1);
if (abs(r - r0*ones(size(r))) > 1e-10)
    error('All points must be at the same radius');
end

theta = atan2(points(:,2), points(:,1));
z = points(:,3);

z = flipud(unique(z));
Nz = length(z);

if (Nz == 1)
    Ntheta = length(theta);
else
    theta1 = theta(1);
    n = 2;
    while (theta(n) ~= theta1)
        n = n + 1;
    end
    Ntheta = n - 1;
end

dtheta = theta(2) - theta(1);
theta = 0:dtheta:(2*pi - dtheta);

eta = reshape(eta0, Ntheta, Nz).';