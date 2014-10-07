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
function [isIn] = isInsideLoop(loop, pnt)

N = size(loop);
if (N(2) ~= 2)
    error('The loop must be an Nx2 matrix of (x,y) points');
end
N = N(1);

outpnt = max(loop) + [1 1];

nCross = 0;

% for debugging
% pla = false;
% if (abs(pnt(1)+0.9)<1e-2 && abs(pnt(2) + 8.7)<1e-2)
%     pla = true;
%     ghljk = 89;
% end

for n = 2:N
    if (segmentsIntersect(loop(n-1,:), loop(n,:), pnt, outpnt))
        nCross = nCross + 1;
%         int = true;
%     else
%         int = false;
    end
%     if (pla)
%         if (int)
%             col = 'k';
%             
%         else
%             col = 'b';
%         end
%         plot([loop(n-1,1), loop(n,1)], [loop(n-1,2), loop(n,2)], col);
%         hold on;
%         plot([pnt(1) outpnt(1)], [pnt(2) outpnt(2)], col);
%     end
end

% for debugging
    R = sqrt(pnt(1)^2 + pnt(2)^2);
    if (R > 15)
        ghlohj = 7;
    end

if (mod(nCross, 2) == 1)
    isIn = true;
else
    isIn = false;
end

