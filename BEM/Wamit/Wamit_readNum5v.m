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
function [v_rad, v_diff, centers] = Wamit_readNum5v(folderpath, runname, bodynames, T, Nbeta, Ndof, g, varargin)
% Reads WAMIT .5v output file
% Returns the radiation velocity (v_rad) and the diffraction velocity 
% (v_diff) at the field points (points)

useSing = 0; % use single percision

if (~isempty(varargin))
    useSing = varargin{1};
end

opts = checkOptions({'bpo'}, varargin);
readBpo = opts(1);

if readBpo
    [centers, ~, Npoints] = Wamit_readBpo(folderpath, bodynames, false);
else
    [centers, ~, ~, Npoints] = Wamit_readPnl(folderpath, bodynames, false);
end

Np = sum(Npoints);

fidx = fopen(fullfile(folderpath, [runname, '.5vx']));
fidy = fopen(fullfile(folderpath, [runname, '.5vy']));
fidz = fopen(fullfile(folderpath, [runname, '.5vz']));
% ignore the header lines
fgetl(fidx);
fgetl(fidy);
fgetl(fidz);

Nper = length(T);


if (useSing)
    v_rad = single(zeros(Nper, Ndof, 3, Np));
    v_diff = single(zeros(Nper, Nbeta, 3, Np));
else
    v_rad = zeros(Nper, Ndof, 3, Np);
    v_diff = zeros(Nper, Nbeta, 3, Np);
end


omega = 2*pi./T;
vdim = g./omega;

for l = 1:Nper
    for n = 1:Np
        fscanf(fidx,'%f',[1 3]);
        fscanf(fidy,'%f',[1 3]);
        fscanf(fidz,'%f',[1 3]);
        for o = 1:Ndof
             % impicitly multiplying by i
             [re_im] = fscanf(fidx,'%f',[1 2]);
             v_rad(l, o, 1, n) = vdim(l)*complex(-re_im(2), re_im(1));
             [re_im] = fscanf(fidy,'%f',[1 2]);
             v_rad(l, o, 2, n) = vdim(l)*complex(-re_im(2), re_im(1));
             [re_im] = fscanf(fidz,'%f',[1 2]);
             v_rad(l, o, 3, n) = vdim(l)*complex(-re_im(2), re_im(1));
        end
    end
    for m = 1:Nbeta
        for o = 1:Np
            fscanf(fidx,'%f',[1 4]);
            fscanf(fidy,'%f',[1 4]);
            fscanf(fidz,'%f',[1 4]);
            % impicitly multiplying by i
            [re_im] =  fscanf(fidx,'%f',[1 2]);
            v_diff(l, m, 1, o) = vdim(l)*complex(-re_im(2), re_im(1));
            [re_im] =  fscanf(fidy,'%f',[1 2]);
            v_diff(l, m, 2, o) = vdim(l)*complex(-re_im(2), re_im(1));
            [re_im] =  fscanf(fidz,'%f',[1 2]);
            v_diff(l, m, 3, o) = vdim(l)*complex(-re_im(2), re_im(1));
        end
    end    
end

fclose(fidx);
fclose(fidy);
fclose(fidz);