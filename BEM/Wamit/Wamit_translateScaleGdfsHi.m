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
function [verts] = Wamit_translateScaleGdfsHi(sourceFolder, gdfIn, destFolder, delta, gdfOut, scale)

if nargin < 5
    gdfOut = gdfIn;
end

if nargin < 6
    scale = 1;
end

if ~iscell(gdfIn)
    gdfIn = {gdfIn};
end

if ~iscell(gdfOut)
    gdfOut = {gdfOut};
end

Ngdf = length(gdfIn);

header = cell(Ngdf, 1);
ulen = zeros(Ngdf, 1);
g = zeros(Ngdf, 1);
isx = zeros(Ngdf, 1);
isy = zeros(Ngdf, 1);
Npat = zeros(Ngdf, 1);
igdf = zeros(Ngdf, 1);

verts = cell(Ngdf, 1);
Ng = cell(Ngdf, 1);
Nk = cell(Ngdf, 1);
knts = cell(Ngdf, 1);

for m = 1:Ngdf
    fid = fopen([sourceFolder '\' gdfIn{m} '.gdf']);
    
    header{m} = fgetl(fid);
    num = textscan(fid,'%f',2);
    num = num{1};
    ulen(m) = num(1);
    g(m) = num(2);

    fgetl(fid);
    num = textscan(fid,'%f',2);
    num = num{1};
    isx(m) = num(1);
    isy(m) = num(2);

    fgetl(fid);
    num = textscan(fid,'%f',2);
    num = num{1};
    Npat(m) = num(1);
    igdf(m) = num(2);
    fgetl(fid);
    
    Ng{m} = zeros(Npat(m), 2);
    Nk{m} = zeros(Npat(m), 2);
    knts{m} = cell(Npat(m), 2);
    verts{m} = cell(Npat(m), 1);
    
    for n = 1:Npat(m)
        num = textscan(fid,'%f',4);
        Nug = num{1}(1);
        Nvg = num{1}(2);
        Kug = num{1}(3);
        Kvg = num{1}(4);
        
        Ng{m}(n,1) = Nug;
        Ng{m}(n,2) = Nvg;
        Nk{m}(n,1) = Kug;
        Nk{m}(n,2) = Kvg;

        Nua = Nug+2*Kug-1;
        Nva = Nvg+2*Kvg-1;
        
        knts{m}(n,1) = textscan(fid,'%f',Nua);
        knts{m}(n,2) = textscan(fid,'%f',Nva);

        Nb = (Nug + Kug - 1)*(Nvg + Kvg - 1);

        num = textscan(fid, '%f', Nb*3);

        verts{m}{n} = reshape(num{1},3,Nb)';
    end
    
    fclose(fid);
end

if ~all(ulen == ulen(1)) || ~all(g == g(1))...
        || ~all(isx == isx(1)) || ~all(isy == isy(1)) 
    warning('Gdf settings are not the same.');
end

for m = 1:Ngdf
    filename = [destFolder '\' gdfOut{m} '.gdf'];
    fid = fopen(filename, 'wt');

    fprintf(fid, header{m});
    fprintf(fid, '\n%9.5f %9.5f\tULEN GRAV\n', ulen(m), g(m));
    fprintf(fid, '%i %i\tISX  ISY\n', isx(m), isy(m));
    fprintf(fid, '%i %i\tNPATCH, IGDEF\n', Npat(m), igdf(m));

    for n = 1:Npat(m)
        
        fprintf(fid, '%i %i\n', Ng{m}(n,1), Ng{m}(n,2));
        fprintf(fid, '%i %i\n', Nk{m}(n,1), Nk{m}(n,2));
        
        for ii = 1:2
            for i = 1:length(knts{m}{n,ii})
                fprintf(fid, ' %8.4f', knts{m}{n,ii}(i));
            end
            fprintf(fid, '\n');
        end
        
        Nv = size(verts{m}{n}, 1);
        for ii = 1:Nv
            for i = 1:3
                fprintf(fid, ' %8.4f', scale*(verts{m}{n}(ii,i)) + delta(i));
            end
            fprintf(fid, '\n');
        end
        fprintf(fid, '\n');
    end

    fclose(fid);
end

end