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
function [err] = matErr(A, As, varargin)

[M1 N1] = size(A);
[M2 N2] = size(As);

if ((M1 ~= M2) || (N1 ~= N2))
    error('matrices must be the same size');
end

if (~isempty(varargin))
   pwr10 = varargin{1};
   A = round2val(A, pwr10);
   As = round2val(As, pwr10);
end

err = abs((A-As)./As);

for m = 1:M1
    for n = 1:N1
        if (isnan(err(m,n)) || isinf(err(m,n)))
            err(m,n) = 0;
        end
    end
end


