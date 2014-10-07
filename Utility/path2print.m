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
function [pathout] = path2print(pathin)
% Takes a normal path string and turns it into a path string that can be
% written with fprintf, sprintf, etc. by adding an extra backslash.  i.e.
% C:\thispath\thisfile.m becomes C:\\thispath\\thisfile.m

folderlist = cell(1,1);

[pathstr name ext] = fileparts(pathin);
folderlist{1} = [name ext];

n = 2;
while (1)
    [pathstr name] = fileparts(pathstr);
    if (~isempty(name))
        folderlist{n} = name;
        n = n + 1;
    else
        break;
    end
end

N = length(folderlist);
pathout = [pathstr '\' folderlist{N}];

for n = N-1:-1:1
    pathout = [pathout '\\' folderlist{n}];
end
