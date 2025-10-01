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
function [centers, areas, norms, Npoints] = Wamit_readPnl(folderpath, bodynames, combineBodies)
% Reads WAMIT .pnl output file
% Returns the point centers, areas, and normals

centers = cell(size(bodynames));
areas = cell(size(bodynames));
norms = cell(size(bodynames));
Npoints = zeros(size(bodynames));

for n = 1:length(bodynames)
    buffer = importdata(fullfile(folderpath, [bodynames{n},'.pnl']));

    raw = buffer.data;

    symCount = length(unique(raw(:,1)));
    if (symCount > 1)
        error('Not set up for symmetry yet...');
    end
    centers{n} = raw(:, 3:5);
    [Npoints(n), ~] = size(centers{n});
    areas{n} = raw(:, 6);
    norms{n} = raw(:, 7:9);
end

if combineBodies
    c1 = centers;
    a1 = areas;
    n1 = norms;
    ind = 0;
    for n = 1:length(bodynames)
        centers((ind+1):(ind+Npoints(n)),:) = c1{n};
        areas((ind+1):(ind+Npoints(n))) = a1{n};
        norms((ind+1):(ind+Npoints(n)),:) = n1{n};
        ind = ind + Npoints(n);
    end
    Npoints = sum(Npoints);
end