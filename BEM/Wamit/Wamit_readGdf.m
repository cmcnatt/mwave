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
function [panGeo] = Wamit_readGdf(folderpath, filename)

fid = fopen(fullfile(folderpath, [filename, '.gdf']));
fgetl(fid);
num = textscan(fid,'%f',2);
num = num{1};
ulen = num(1);
g = num(2);

fgetl(fid);
num = textscan(fid,'%f',2);
num = num{1};
isx = num(1);
isy = num(2);

fgetl(fid);
num = textscan(fid,'%f',1);
Npan = num{1};

pans(Npan) = Panel;

fgetl(fid);
for n = 1:Npan
    num = textscan(fid,'%f',12);
    verts = reshape(num{1},3,4)';
    
    pans(n) = Panel(verts);
    fgetl(fid);
end

if (~isx && ~isy)
    panGeo = PanelGeo(pans);
elseif (isx && ~isy)
    panGeo = PanelGeo(pans, 'xsym');
elseif (~isx && isy)
    panGeo = PanelGeo(pans, 'ysym');
elseif (isx && isy)
    panGeo = PanelGeo(pans, 'xsym', 'ysym');
else
    error('Symmetry values not recognized');
end

fclose(fid);