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
function [panGeo] = Wamit_readGdfHi(folderpath, filename)

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
num = textscan(fid,'%f',2);
num = num{1};
Npan = num(1);
igdef = num(2);
fgetl(fid);

pans(Npan) = Panel;

if igdef == 0
    for n = 1:Npan
        num = textscan(fid,'%f',12);
        verts = reshape(num{1},3,4)';

        pans(n) = Panel(verts);
        fgetl(fid);
    end
elseif igdef == 1
    verts = cell(Npan, 1);
    Knts = cell(Npan, 2);

    for n = 1:Npan
        num = textscan(fid,'%f',4);
        Nug = num{1}(1);
        Nvg = num{1}(2);
        Kug = num{1}(3);
        Kvg = num{1}(4);

        Nua = Nug+2*Kug-1;
        Nva = Nvg+2*Kvg-1;
        Knts(n,1) = textscan(fid,'%f',Nua);
        Knts(n,2) = textscan(fid,'%f',Nva);

        Nb = (Nug + Kug - 1)*(Nvg + Kvg - 1);

        num = textscan(fid, '%f', Nb*3);

        vertsn = reshape(num{1},3,Nb)'; 
        verts{n} = vertsn;          % verts in WAMIT order  
        vertsn4 = vertsn(4, :);     % need to swap points 3 and 4 to make panel
        vertsn(4, :) = vertsn(3, :);
        vertsn(3, :) = vertsn4;

        pans(n) = Panel(vertsn);
    end
end

panGeo = PanelGeo(pans);

fclose(fid);

end