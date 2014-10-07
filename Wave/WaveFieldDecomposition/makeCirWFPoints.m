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
function [points] = makeCirWFPoints(r, Ntheta, z)
% Creates the circular points needed for WAMIT in order to conduct a
% cylindrical wave field decomposition.
%
%   - r: The radius at which to compute the wave decomposition
%   - Ntheta: The number of points in the circular direction. The points 
%   will be evenly spaced around the circel. The more points the more 
%   accurate the wave field decomposition.  Using a number that is a power 
%   of 2 will speed up the FFT in the wave field decomposition.
%   - The points in the z-direction. It is recommended to use cosine
%   spacing near the z = 0 free-surface to improve resolution. For example,
%       
%       z = -h*(1-cos(0:pi/2/nZ:pi/2));
%       z = round(z*10^6)/10^6;
%
%   The round is used to ensure that the first point is exactly at z = 0.
%   i.e. to get rid of rounding errors.

Nz = length(z);
Nr = length(r);
points = zeros(Nr*Ntheta*Nz, 3);

delTheta = 2*pi/Ntheta;
theta = 0:delTheta:(2*pi - delTheta);
ctheta = cos(theta);
stheta = sin(theta);

for m = 1:Nr   
    for n = 1:Nz
        start = (m-1)*Nz*Ntheta + (n-1)*Ntheta + 1;
        range = (start:start + Ntheta - 1);
        
        points(range, 1) = r(m)*ctheta;
        points(range, 2) = r(m)*stheta;
        points(range, 3) = z(n)*ones(Ntheta, 1);
    end
end