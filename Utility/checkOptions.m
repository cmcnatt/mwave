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
function [opts, args] = checkOptions(posOpts, inputs)
% Utility function to check varargin options
% Assumes the varargin comes in the form: 'OptionName', arg1, arg2, etc.,
% where arg1, arg2 are arguments associated with the option called
% 'OptionName'
%
% INPUTS:
% posOpts - a cell array of all the possible options that the 
% function/method could take. 
% inputs - the actual inputs, should just be: varargin
%
% OUTPUTS:
% opts - a boolean array of the same length as the number of possible
% options, which indicated whether that possible option was selected.
% args - a cell array of the arguements, where the first cell index
% corresponds to the option and the second index corresponds to index of
% the argument.
%
% EXAMPLE: 
% [opts, args] = checkOptions({{'Beef', 2}, {'Wine, 4}, 'Cheese'}, varargin)
% where in this case, varargin is 'Cheese', 'Beef', true, 3
% the outputs are:
% opts = [1 0 1] to indicate that the 'Beef' and 'Cheese' options had been
% selected (the first and third in the posOpts list)
% args = {{true, 3}, -1, -1} where args{1} contains the cell array of the 
% 'Beef' arguements, and args{2} = -1 because the 'Wine' option was not
% selected and args{3} = -1 because 'Cheese' didn't have any arguments.

N = length(inputs);
M = length(posOpts);

opts = false*ones(1, M);
args = cell(1, M);

for m = 1:M
    posOptm = posOpts{m};
    nargOptm = 0;
    if (iscell(posOptm))
        if (length(posOptm) > 1)
            nargOptm = posOptm{2};
        end
        posOptm = posOptm{1};
    end
    
    % default value
    args{m} = NaN;

    for n = 1:N
        inputn = inputs{n};
        if (strcmpi(inputn, posOptm))
            opts(m) = true;
            if (nargOptm > 0)
                if (nargOptm == 1)
                    args{m} = inputs{n+1};
                else
                    argsm = cell(1, nargOptm);
                    for l = 1:nargOptm
                        argsm{l} = inputs{n+l};
                    end
                    args{m} = argsm;
                end
            end
        end
    end
end

