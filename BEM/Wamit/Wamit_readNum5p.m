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
function [p_rad, p_diff, centers] = Wamit_readNum5p(folderpath, runname, bodynames, Nper, Nbeta, Ndof, rho, g, cgs, varargin)
% Reads WAMIT .5p output file
% Returns the radiation pressure (p_rad) and the diffraction pressure 
% (p_diff) at the field points (points)

useSing = 0; % use single percision

if (~isempty(varargin))
    useSing = varargin{1};
end

opts = checkOptions({{'bpo'} {'ihi'}}, varargin);
readBpo = opts(1);
ihi = opts(2);

if readBpo
    [centers, ~, Npoints, noConv] = Wamit_readBpo(folderpath, bodynames, false, cgs, ihi);
    Np = sum(Npoints);
    noConv = cell2mat(noConv.');
else
    [centers, ~, ~, Npoints] = Wamit_readPnl(folderpath, bodynames, false);
    Np = sum(Npoints);
    noConv = zeros(Np, 1);
end

if (useSing)
    p_rad = single(zeros(Nper, Ndof, Np));
    p_diff = single(zeros(Nper, Nbeta, Np));
else
    p_rad = zeros(Nper, Ndof, Np);
    p_diff = zeros(Nper, Nbeta, Np);
end

pdim = rho*g;
fid = fopen(fullfile(folderpath, [runname, '.5p']));
% ignore the header line
fgetl(fid);

for l = 1:Nper
    for n = 1:Np
        buffer = fscanf(fid,'%f',[1 3]); 
        for o = 1:Ndof
            [re_im] = fscanf(fid,'%f',[1 2]);
            if noConv(n)
                p_rad(l, o, n) = 0;
            else
                p_rad(l, o, n) = pdim*complex(re_im(1), re_im(2));
            end
        end
    end
    for m = 1:Nbeta
        for o = 1:Np
            buffer = fscanf(fid,'%f',[1 4]);
            [re_im] =  fscanf(fid,'%f',[1 2]);
             if noConv(o)
                 p_diff(l, m, o) = 0;
             else
                 p_diff(l, m, o) = pdim*complex(re_im(1), re_im(2));
             end
        end
    end    
end

fclose(fid);