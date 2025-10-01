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
function [] = Wamit_writeBpi(fileLoc, name, points)

if isempty(points)
    Np = 0;
else
    [Np, M] = size(points);
    if 3 ~= M
        error('points must be an Np x 3 matrix');
    end
end
    
fileName = fullfile(fileLoc,[name, '.bpi']);
fileID = fopen(fileName, 'wt');

fprintf(fileID, ['Model ' name ', created: ' date '\n']);
fprintf(fileID, '%i\n', Np);
for n = 1:Np
    fprintf(fileID, '\t%8.4f\t%8.4f\t%8.4f\n', points(n,1), points(n,2), points(n,3));
end

fclose(fileID);
end