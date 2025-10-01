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
function [C] = Wamit_readHst(folderpath, runname, rho, g)

file_data = importdata(fullfile(folderpath, [runname, '.hst']));
data = file_data.data;

dof = max(data(:,1));

C = zeros(dof, dof);

for n = 1:dof
    for m = 1:dof
        C(n,m) = data((n-1)*dof + m, 3);
    end
end

C = rho*g*C;