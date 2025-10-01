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
function [v_rad, v_diff, points] = Wamit_readNum6v(folderpath, runname, T, Nbeta, Ndof, g)
% Reads WAMIT .7 output file
% Returns the radiation velocity (v_rad) and the diffraction velocity 
% (v_diff) at the field points (points)

points = Wamit_readFpt(folderpath, runname);

[Npoints, buffer] = size(points);

fidx = fopen(fullfile(folderpath, [runname, '.6vx']));
fidy = fopen(fullfile(folderpath, [runname, '.6vy']));
fidz = fopen(fullfile(folderpath, [runname, '.6vz']));
% ignore the header lines
fgetl(fidx);
fgetl(fidy);
fgetl(fidz);

Nper = length(T);

v_rad = zeros(Nper, Ndof, 3, Npoints);
v_diff = zeros(Nper, Nbeta, 3, Npoints);

omega = 2*pi./T;
vdim = g./omega;

for l = 1:Nper
    for m = 1:Npoints
        buffer = fscanf(fidx,'%f',[1 2]);
        buffer = fscanf(fidy,'%f',[1 2]);
        buffer = fscanf(fidz,'%f',[1 2]);
        for n = 1:Ndof
             % impicitly multiplying by i
             [re_im] = fscanf(fidx,'%f',[1 2]);
             v_rad(l, n, 1, m) = vdim(l)*complex(-re_im(2), re_im(1));
             [re_im] = fscanf(fidy,'%f',[1 2]);
             v_rad(l, n, 2, m) = vdim(l)*complex(-re_im(2), re_im(1));
             [re_im] = fscanf(fidz,'%f',[1 2]);
             v_rad(l, n, 3, m) = vdim(l)*complex(-re_im(2), re_im(1));
        end
    end
    for m = 1:Nbeta
        for n = 1:Npoints
            buffer = fscanf(fidx,'%f',[1 3]);
            buffer = fscanf(fidy,'%f',[1 3]);
            buffer = fscanf(fidz,'%f',[1 3]);
            % impicitly multiplying by i
            [re_im] =  fscanf(fidx,'%f',[1 2]);
            v_diff(l, m, 1, n) = vdim(l)*complex(-re_im(2), re_im(1));
            [re_im] =  fscanf(fidy,'%f',[1 2]);
            v_diff(l, m, 2, n) = vdim(l)*complex(-re_im(2), re_im(1));
            [re_im] =  fscanf(fidz,'%f',[1 2]);
            v_diff(l, m, 3, n) = vdim(l)*complex(-re_im(2), re_im(1));
        end
    end    
end

fclose(fidx);
fclose(fidy);
fclose(fidz);