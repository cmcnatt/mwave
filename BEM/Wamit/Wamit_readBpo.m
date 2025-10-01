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
function [centers, norms, Npoints, noConv] = Wamit_readBpo(folderpath, bodynames, combineBodies, cgs, ilowHi)
% Reads WAMIT .bpo output file
% Returns the point centers, and normals

if nargin < 4
    ilowHi = 1;
end

centers = cell(size(bodynames));
norms = cell(size(bodynames));
Npoints = zeros(size(bodynames));
noConv = cell(size(bodynames));

totPanCnt = 0;
for n = 1:length(bodynames)
    name = fullfile(folderpath,[bodynames{n},'.bpo']);
    if exist(name, 'file')
        buffer = importdata(name, ' ', 2);

        raw = buffer.data;

        if ilowHi
            centers{n} = raw(:, 8:10);
            [Npoints(n), ~] = size(centers{n});
            norms{n} = raw(:, 11:13);
            it = raw(:,7);
            noConv{n} = it > 16;
        else
            bodyGeo = Wamit_readGdf(folderpath, bodynames{n});
            bodyGeo.Translate(cgs{n});
            ipan = raw(:,2);
            ipan = ipan - totPanCnt;
            Npoints(n) = length(ipan);

            centers{n} = bodyGeo.Centroids(ipan,:);
            norms{n} = bodyGeo.Normals(ipan,:);
            noConv{n} = zeros(Npoints(n), 1);
        end
    else
        if ~ilowHi
            bodyGeo = Wamit_readGdf(folderpath, bodynames{n});
        end
    end
    if ~ilowHi
        totPanCnt = totPanCnt + bodyGeo.Count;
    end
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