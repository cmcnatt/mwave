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
function [A, B, T, Modes, Ainf, A0] = Wamit_read1(folderpath, runname, rho, varargin)
% reads WAMIT .1 output file
% returns the added mass (A), damping (B), the periods (T), the Modes 
% (Modes), infinite frequency, zero period added mass (Ainf),
% zero frequency, infinite period added mass (A0);

% Check if spikes should also be removed from coefficients (based on
% user-specified frequency ranges, e.g. a spike due to a trapped mode resulting from the mesh design).
[opts, args] = checkOptions({{'removeSpikes',1}}, varargin);

if opts(1)
    spikeFreqs = args{1};
    if ~iscell(spikeFreqs)
        error('SpikeFreqs must be a cell array.')
    end
    if size(spikeFreqs,1) ~= 2
        error('SpikeFreqs must contain two rows - the first for in-plane modes, the second for out-of-plane modes.')
    end
    for i = 1:size(spikeFreqs,1)
        for j = 1:size(spikeFreqs,2)
            if ~(size(spikeFreqs{i,j},1)==1 && size(spikeFreqs{i,j},2)==2)
                error('Each cell in cell array must be a 1x2 vector containing the lower and upper frequency bounds.')
            end
        end
    end
end

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

if opts(1)
    %%  Use interpolation to smooth out spikes (if desired)
    % Find indices that correspond to the frequency ranges provided
    warning('Remember that the spike removal option is only set up to work with a two body WAMIT run (planar or 6DoF).')
    w = 2*pi./T2;
    for i = 1:2 % Row 1 contains in-plane freq ranges, row 2 out-of-plane freq ranges.
        for j = 1:size(spikeFreqs,2) % No. of spikes to interpolate out
            for k = 1:2 % upper and lower bounds
                [~,interpInds{i,j}(k)] = min(abs(w - spikeFreqs{i,j}(k)));
            end
            if interpInds{i,j}(1) == interpInds{i,j}(2)
                warning('Upper bound index == lower bound index: This may suggest user did not run WAMIT over enough frequencies, or that the spike frequencies are too close together.')
            end
        end
    end
    if length(Modes)==12
        DOFinds{1,1} = [1,3,5,7,9,11]; % in-plane DoFs
        DOFinds{2,1} = [2,4,6,8,10,12]; % out-of-plane DoFs
    elseif length(Modes)==6
        DOFinds{1,1} = [1:6]; % All are in-plane DoFs
    else
        error('No. of modes in WAMIT run is not compatible with spike removal option.')
    end
    for i = 1:round(length(Modes)/6) % i=1: remove spikes from in-plane DoFs, i=2: from out-of-plane DoFs
        for j = 1:size(spikeFreqs,2)
            interp_lb = interpInds{i,j}(1);
            interp_ub = interpInds{i,j}(2);
            for k = interp_lb:interp_ub
                if interp_ub == length(w)
                    B(k,DOFinds{i,1},DOFinds{i,1}) = B(interp_lb-1,DOFinds{i,1},DOFinds{i,1});
                    A(k,DOFinds{i,1},DOFinds{i,1}) = A(interp_lb-1,DOFinds{i,1},DOFinds{i,1});
                else
                    B(k,DOFinds{i,1},DOFinds{i,1}) = B(interp_lb-1,DOFinds{i,1},DOFinds{i,1}) + ...
                        (B(interp_ub+1,DOFinds{i,1},DOFinds{i,1})-B(interp_lb-1,DOFinds{i,1},DOFinds{i,1})).*(k-(interp_lb-1))./((interp_ub+1)-(interp_lb-1));
                    A(k,DOFinds{i,1},DOFinds{i,1}) = A(interp_lb-1,DOFinds{i,1},DOFinds{i,1}) + ...
                        (A(interp_ub+1,DOFinds{i,1},DOFinds{i,1})-A(interp_lb-1,DOFinds{i,1},DOFinds{i,1})).*(k-(interp_lb-1))./((interp_ub+1)-(interp_lb-1));
                end
                if interp_lb == 1
                    error('Code not currently set up to use lowest frequency as part of interpolation range.\n%s',...
                        'Spikes are not likely to occur here anyway if choose frequency range prudently.')
                end
            end
        end
    end
    
    if length(Modes)==12 %i.e. each body has 6DoF
        % For couplings between in-plane and out-of-plane degrees of freedom
        for i = 1:2 % i=1: remove spikes from couplings to in-plane from out-of-plane DoFs, i=2: from couplings to out-of-plane from in-plane.
            for j = 1:size(spikeFreqs,2)
                interp_lb = interpInds{i,j}(1);
                interp_ub = interpInds{i,j}(2);
                for k = interp_lb:interp_ub
                    if interp_ub == length(w)
                        B(k,DOFinds{i,1},DOFinds{3-i,1}) = B(interp_lb-1,DOFinds{i,1},DOFinds{3-i,1});
                        A(k,DOFinds{i,1},DOFinds{3-i,1}) = A(interp_lb-1,DOFinds{i,1},DOFinds{3-i,1});
                    else
                        B(k,DOFinds{i,1},DOFinds{3-i,1}) = B(interp_lb-1,DOFinds{i,1},DOFinds{3-i,1}) + ...
                            (B(interp_ub+1,DOFinds{i,1},DOFinds{3-i,1})-B(interp_lb-1,DOFinds{i,1},DOFinds{3-i,1})).*(k-(interp_lb-1))./((interp_ub+1)-(interp_lb-1));
                        A(k,DOFinds{i,1},DOFinds{3-i,1}) = A(interp_lb-1,DOFinds{i,1},DOFinds{3-i,1}) + ...
                            (A(interp_ub+1,DOFinds{i,1},DOFinds{3-i,1})-A(interp_lb-1,DOFinds{i,1},DOFinds{3-i,1})).*(k-(interp_lb-1))./((interp_ub+1)-(interp_lb-1));
                    end
                end
            end
        end
    end
    
end

T = T2;










