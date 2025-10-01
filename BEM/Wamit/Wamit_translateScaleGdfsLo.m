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
function [] = Wamit_translateScaleGdfsLo(sourceFolder, gdfIn, destFolder, delta, gdfOut, scale)

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

Npan = zeros(Ngdf, 1);

ulen = zeros(Ngdf, 1);
g = zeros(Ngdf, 1);
isx = zeros(Ngdf, 1);
isy = zeros(Ngdf, 1);

verts = cell(Ngdf, 1);

for n = 1:Ngdf
    fid = fopen(fullfile(sourceFolder, [gdfIn{n}, '.gdf']));
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

for m = 1:Ngdf
    filename = [destFolder '\' gdfOut{m} '.gdf'];
    fid = fopen(filename, 'wt');

    fprintf(fid, ['Model ' gdfOut{m} ', created: ' date '\n']);
    fprintf(fid, '%9.5f %9.5f\n', ulen(m), g(m));
    fprintf(fid, '%i %i \n', isx(m), isy(m));
    fprintf(fid, '%i\n', Npan(m));

    i = 1;
    for n = 1:length(verts{m})
        fprintf(fid, '\t%9.5f', scale*verts{m}(n) + delta(i));
        i = i + 1;
        if i > 3
            i = 1;
            fprintf(fid, '\n');
        end
    end

    fclose(fid);
end



