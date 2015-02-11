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
function [eta, points] = nemoh_readFieldPoints(fullPath, isFS)

fid = fopen(fullPath);

fgetl(fid); % ignore the fist header line
lin = fgetl(fid);

npts = str2double(lin(8:15));

if (isFS)
    pntsInd = 2;
else
    pntsInd = 3;
end

points = zeros(npts, 3);
eta = zeros(npts, 1);

for n = 1:npts
    lin = fscanf(fid,'%f', pntsInd+4);
    points(n,1:pntsInd) = lin(1:pntsInd);
    % Note that the complex conjuage is taken here because Nemoh
    % computes for a time depedence of exp(-i*omega*t), while
    % mwave assumes a time depedent of exp(i*omega*t);
    eta(n) = lin(pntsInd+1)*exp(-1i*lin(pntsInd+2));
end

fclose(fid);
end