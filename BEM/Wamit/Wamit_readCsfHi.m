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
    C. McNatt, A. Cotten
%}
function [panGeoCsf] = Wamit_readCsfHi(folderpath, filename)

fid = fopen(fullfile(folderpath, [filename,'.csf']);
fgetl(fid);
num = textscan(fid,'%f',1);
num = num{1};
iLowHiCsf = num(1);

fgetl(fid);
num = textscan(fid,'%f',2);
num = num{1};
isxcsf = num(1);
isycsf = num(2);

fgetl(fid);
num = textscan(fid,'%f',3);
num = num{1};
Npatcsf = num(1);
icdef = num(2);
pszcsf = num(3);
fgetl(fid);

pans(Npatcsf) = Panel;

if icdef == 0
    for n = 1:Npatcsf
        num = textscan(fid,'%f',12);
        verts = reshape(num{1},3,4)';

        pans(n) = Panel(verts);
        fgetl(fid);
    end
elseif icdef == 1
    verts = cell(Npatcsf, 1);
    Knts = cell(Npatcsf, 2);

    for n = 1:Npatcsf
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

panGeoCsf = PanelGeo(pans);

fclose(fid);

end