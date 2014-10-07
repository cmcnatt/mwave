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
function [] = addMwaveHeader(fileLoc, fileName, varargin)
% Add the mwave header (found in the script file: mwaveHeader) to a 
% specific file. If the file already has the header (it check the first 
% line), it does not add it.
%
% The optional input argument, 'OnlyM', ensures that the header is added to
% only input files with the .m extention

opts = checkOptions({'OnlyM'}, varargin);

fileOk = true;
if (opts(1))
    idot = find(fileName == '.');
    if (strcmp(fileName(idot:end), '.m'))
        fileOk = true;
    else
        fileOk = false;
    end
end

if (fileOk)
    fidhead = fopen([mwavePath '\Header\mwaveHeader.m'], 'rt');
    header = [];
    n = 1;
    while true
        if (n == 1)
            head1 = fgets(fidhead);
            thisline = head1;
        else
            thisline = fgets(fidhead);
        end
        if ~ischar(thisline); break; end   %end of file
        header = [header thisline];
        n = n + 1;
    end
    fclose(fidhead);

    fid = fopen([fileLoc '\' fileName], 'rt');
    line1 = fgets(fid);
    fclose(fid);

    if (~strcmp(head1, line1))
        addHeader(fileLoc, fileName, header);
    end

    if (exist([fileLoc '\' fileName '~'], 'file'))
        delete([fileLoc '\' fileName '~']);
    end
end


