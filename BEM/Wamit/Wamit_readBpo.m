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
function [centers, norms, Npoints, noConv] = Wamit_readBpo(folderpath, bodynames, combineBodies)
% Reads WAMIT .bpo output file
% Returns the point centers, and normals

centers = cell(size(bodynames));
norms = cell(size(bodynames));
Npoints = zeros(size(bodynames));
noConv = cell(size(bodynames));

for n = 1:length(bodynames)
    buffer = importdata([folderpath '\' bodynames{n} '.bpo'], ' ', 2);

    raw = buffer.data;

    centers{n} = raw(:, 8:10);
    [Npoints(n), ~] = size(centers{n});
    norms{n} = raw(:, 11:13);
    it = raw(:,7);
    noConv{n} = it > 16;
end

if combineBodies
    c1 = centers;
    n1 = norms;
    nc = noConv;
    ind = 0;
    for n = 1:length(bodynames)
        centers((ind+1):(ind+Npoints(n)),:) = c1{n};
        norms((ind+1):(ind+Npoints(n)),:) = n1{n};
        noConv((ind+1):(ind+Npoints(n))) = nc{n};
        ind = ind + Npoints(n);
    end
    Npoints = sum(Npoints);
end