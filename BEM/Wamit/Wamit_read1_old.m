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
function [A, B, T, Modes] = Wamit_read1(folderpath, runname, rho)
% reads WAMIT .1 output file
% returns the added mass (A), damping (B), the periods (T), and the Modes 
% (Modes)

% Read in the file and ignore the header line.
file_data = importdata(fullfile(folderpath, [runname,'.1']));
data = file_data.data;

% First, find out how many periods there are and make a vector of them.
T = unique(data(:,1));

if (T(1) ~= data(1,1));
    T = flipud(T);
end

for n = 1:length(T)
    if (T(n) < 0)
        T(n) = Inf;
    end
end

% for now, remove zero and infinite frequency results

% how many modes are there.  make a vector of the indecies.
Modes1 = unique(data(:,2));
m = 1;
for n = 1:length(Modes1)
    if (~isnan(Modes1(n)))
        Modes(m) = Modes1(n);
        m = m + 1;
    end
end

nT = length(T);
nM = length(Modes);

A = zeros(nT, nM, nM);
B = zeros(nT, nM, nM);

% for i = 1:length(T)
%     II = find(data(:,1)==T(i));    
%     for j = 1:length(II)
%         A(i,data(II(j),2),data(II(j),3)) = data(II(j),4);
%         B(i,data(II(j),2),data(II(j),3)) = data(II(j),5);
%     end
% end

Omega = 2*pi./T;

lastM = Modes(end);
modei = zeros(lastM,1);
for n = 1:nM
    modei(Modes(n)) = n;
end


N = size(data, 1);
lt = data(1,1);
i = 1;
for n = 1:N  
    t = data(n,1);
    if (abs(t-lt) > 0.0001)
        i = i+1;
    end
    j = modei(data(n, 2));
    k = modei(data(n, 3));
    A(i,j,k) = rho*data(n, 4);
    B(i,j,k) = rho*Omega(i)*data(n, 5);
    lt = t;
end










