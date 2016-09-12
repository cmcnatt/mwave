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
function [A, B, T, Modes, Ainf, A0] = Wamit_read1(folderpath, runname, rho)
% reads WAMIT .1 output file
% returns the added mass (A), damping (B), the periods (T), the Modes 
% (Modes), infinite frequency, zero period added mass (Ainf),
% zero frequency, infinite period added mass (A0);

% Read in the file and ignore the header line.
fid = fopen([folderpath '\' runname '.1']);
head = fgetl(fid);

n = 0;
while ~feof(fid)
    n = n + 1;
    vals = str2num(fgetl(fid));
    modes(n, 1:2) = vals(2:3);
    
    t(n) = vals(1);
    a(n) = vals(4);
    
    if 0 == t(n)
        b(n) = 0;
    elseif 0 > t(n)
        t(n) = Inf;
        b(n) = 0;
    else
        b(n) = vals(5);
    end
end

fclose(fid);

N = n;

T = unique(t, 'stable');
Modes = unique(modes);

nT = length(T);
nM = length(Modes);

A = zeros(nT, nM, nM);
B = zeros(nT, nM, nM);

Omega = 2*pi./t;

lastM = Modes(end);
modei = zeros(lastM,1);
for n = 1:nM
    modei(Modes(n)) = n;
end


lt = t(1);
i = 1;
for n = 1:N  
    if (abs(t(n)-lt) > 0.0001)
        i = i+1;
    end
    j = modei(modes(n,1));
    k = modei(modes(n,2));
    A(i,j,k) = rho*a(n);
    B(i,j,k) = rho*Omega(n)*b(n);
    lt = t(n);
end

Ainf = [];
A0 = [];

m = 0;
for n = 1:nT
    if (Inf == T(n)) || (0 == T(n))
        if Inf == T(n)
            % A0 is the zero frequency, infinite period added mass
            A0 = squeeze(A(m+1,:,:));
        else
            Ainf = squeeze(A(m+1,:,:));
        end
        Alow = A(1:m,:,:);
        Ahigh = A(m+2:end,:,:);
        A = cat(1, Alow, Ahigh);

        Blow = B(1:m,:,:);
        Bhigh = B(m+2:end,:,:);
        B = cat(1, Blow, Bhigh);            
    else
        m = m + 1;
        T2(m) = T(n);
    end
end

T = T2;










