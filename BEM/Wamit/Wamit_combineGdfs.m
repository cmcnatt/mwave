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
function [Npan] = Wamit_combineGdfs(sourceFolder, gdfIn, destFolder, gdfOut)

Ngdf = length(gdfIn);

Npan = zeros(Ngdf, 1);

ulen = zeros(Ngdf, 1);
g = zeros(Ngdf, 1);
isx = zeros(Ngdf, 1);
isy = zeros(Ngdf, 1);

verts = cell(Ngdf, 1);

for n = 1:Ngdf
    fid = fopen(fullfile(sourceFolder,[gdfIn{n},'.gdf']));
    fgetl(fid);
    num = textscan(fid,'%f',2);
    num = num{1};
    ulen(n) = num(1);
    g(n) = num(2);

    fgetl(fid);
    num = textscan(fid,'%f',2);
    num = num{1};
    isx(n) = num(1);
    isy(n) = num(2);

    fgetl(fid);
    num = textscan(fid,'%f',1);
    Npan(n) = num{1};
    
    verts(n) = textscan(fid,'%f',12*Npan(n));
    fclose(fid);
end

if ~all(ulen == ulen(1)) || ~all(g == g(1))...
        || ~all(isx == isx(1)) || ~all(isy == isy(1)) 
    warning('Gdf settings are not the same.');
end

filename = fullfile(destFolder, [gdfOut, '.gdf']);
fid = fopen(filename, 'wt');
                        
fprintf(fid, ['Model ' gdfOut ', created: ' date '\n']);
fprintf(fid, '%9.5f %9.5f\n', ulen(1), g(1));
fprintf(fid, '%i %i \n', isx(1), isy(1));
fprintf(fid, '%i\n', sum(Npan));

for m = 1:Ngdf
    i = 1;
    for n = 1:length(verts{m})
        fprintf(fid, '\t%9.5f', verts{m}(n));
        i = i + 1;
        if i > 3
            i = 1;
            fprintf(fid, '\n');
        end
    end
end
               
fclose(fid);



