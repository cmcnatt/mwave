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
function [geo] = nemoh_readMesh(fullPath)

fid = fopen(fullPath, 'r');

lin = fscanf(fid, '%i',2);
nPnts = lin(1);
nPans = lin(2);

pnts = zeros(nPnts, 3);

for n = 1:nPnts
    pnts(n,:) = fscanf(fid, '%f', 3);
end

pans(nPans,1) = Panel;

for m = 1:nPans
    vertinds = fscanf(fid, '%i', 4);
    verts = zeros(4,3);
    for n = 1:4
        verts(n,:) = pnts(vertinds(n),:);
    end
     pans(m) = Panel(verts);
end

fclose(fid);

geo = PanelGeo(pans);

end