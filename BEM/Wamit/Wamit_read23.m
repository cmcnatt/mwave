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
function [Forces, T, Beta, Modes] = Wamit_read23(folderpath, runname, rho, g, varargin)
% reads WAMIT .2 output file
% returns the diffraction forces (Forces), the periods (T), the headings 
% (Beta), and the Modes (Modes)

% Check if spikes should also be removed from coefficients (based on
% user-specified frequency ranges, e.g. a spike due to a trapped mode resulting from the mesh design).
[opts, args] = checkOptions({{'fk'},{'sc'},{'removeSpikes',1}}, varargin);

if opts(1) && opts(2)
    warning('User should not select for both ''fk'' and ''sc'' excitation forces to be read from WAMIT files by Wamit_read23 function.\n%s','Nonetheless, in this case, scattering forces will be read.')
    ext = 'sc';
elseif opts(1) % 'fk' argument
    ext = 'fk';
elseif opts(2) % 'sc' argument
    ext = 'sc';
else
    % If no option is selected, WAMIT just outputs a .2 or .3 file, so
    % don't need any further letters for file extension.
    ext = '';
end

if opts(3)
    spikeFreqs = args{3};
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

% read in the file and ignore the header line.
try
    file_data = importdata([folderpath '/' runname '.2' ext]);
catch
    file_data = importdata([folderpath '/' runname '.3' ext]);
end
data = file_data.data;

% First, find out how many periods there are and make a vector of them.
T = unique(data(:,1));

% I want the frequencies to go in ascending order, not the periods, so I
% will reverse the order of the periods.
if (T(1) ~= data(1,1));
    T = flipud(T);
end
%%% BUT it doesn't appear that this T vector is even used anymore, so
%%% ignore the above comment about reversing the periods.

% Second, find out how many direction there are and make a vector of them.
Beta = unique(data(:,2));

% how many modes are there.  make a vector of the indecies.
Modes = unique(data(:,3));

nb = length(Beta);
nt = length(T);
dof = length(Modes);

% make and empty matrix to read in the excitation forces
Forces = zeros(nt,nb,dof);

for i = 1:nt
    for j = 1:nb
        for k = 1:dof
            [re_im] = data((i-1)*nb*dof + (j-1)*dof + k, 6:7);
            Forces(i, j, k) = rho*g*complex(re_im(1), re_im(2));
        end
    end
end

if opts(3)
    %%  Use interpolation to smooth out spikes (if desired)
    % Find indices that correspond to the frequency ranges provided
    warning('Remember that the spike removal option is only set up to work with a two body WAMIT run (planar or 6DoF).')
    w = 2*pi./T;
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
                    Forces(k,:,DOFinds{i,1}) = Forces(interp_lb-1,:,DOFinds{i,1});
                else
                    Forces(k,:,DOFinds{i,1}) = Forces(interp_lb-1,:,DOFinds{i,1}) + ...
                        (Forces(interp_ub+1,:,DOFinds{i,1})-Forces(interp_lb-1,:,DOFinds{i,1})).*(k-(interp_lb-1))./((interp_ub+1)-(interp_lb-1));
                end
                if interp_lb == 1
                    error('Code not currently set up to use lowest frequency as part of interpolation range.\n%s',...
                        'Spikes are not likely to occur here anyway if choose frequency range prudently.')
                end
            end
        end
    end
    
end


