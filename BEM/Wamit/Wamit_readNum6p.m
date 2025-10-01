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
function [p_rad, p_diff, points] = Wamit_readNum6p(folderpath, runname, Nper, Nbeta, Ndof, rho, g, varargin)
% Reads WAMIT .6 output file
% Returns the radiation pressure (p_rad) and the diffraction pressure 
% (p_diff) at the field points (points)

useSing = 0; % use single percision

if (~isempty(varargin))
    useSing = varargin{1};
end

points = Wamit_readFpt(folderpath, runname);

[Npoints buffer] = size(points);


fid = fopen(fullfile(folderpath, [runname, '.6p']));
% ignore the header line
fgetl(fid);

if (useSing)
    p_rad = single(zeros(Nper, Ndof, Npoints));
    p_diff = single(zeros(Nper, Nbeta, Npoints));
else
    p_rad = zeros(Nper, Ndof, Npoints);
    p_diff = zeros(Nper, Nbeta, Npoints);
end

pdim = rho*g;

for l = 1:Nper
    for m = 1:Npoints
        if Ndof > 0
            buffer = fscanf(fid,'%f',[1 2]); 
            for n = 1:Ndof
                [re_im] = fscanf(fid,'%f',[1 2]);
                p_rad(l, n, m) = pdim*complex(re_im(1), re_im(2));
            end
        end
    end
    for m = 1:Nbeta
        for n = 1:Npoints
            buffer = fscanf(fid,'%f',[1 3]);
            [re_im] =  fscanf(fid,'%f',[1 2]);
            p_diff(l, m, n) = pdim*complex(re_im(1), re_im(2));
        end
    end    
end

fclose(fid);