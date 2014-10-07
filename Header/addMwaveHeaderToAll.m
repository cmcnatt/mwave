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
function [] = addMwaveHeaderToAll(fileLoc)
% Adds the mwave header (found in the script file: mwaveHeader) to all the 
% files in the directory and all subdirectories in the directory. It will 
% not add the header to files that already have the header.

curdir = cd(fileLoc);

filesFolds = dir;

N = length(filesFolds);

for n = 1:N
    if (filesFolds(n).isdir)
        if ~(strcmp('.', filesFolds(n).name) || (strcmp('..', filesFolds(n).name)))
            addMwaveHeaderToAll([fileLoc '\' filesFolds(n).name]);
        end
    else
        addMwaveHeader(fileLoc, filesFolds(n).name, 'OnlyM');
    end
end

cd(curdir);