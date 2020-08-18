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
function [verts] = Wamit_translateScaleCsfsHi(sourceFolder, csfIn, destFolder, delta, csfOut, scale)
% This function is the CSF file equivalent of the GDF version.



if nargin < 5
    csfOut = csfIn;
end

if nargin < 6
    scale = 1;
end

if ~iscell(csfIn)
    csfIn = {csfIn};
end

if ~iscell(csfOut)
    csfOut = {csfOut};
end

Ngdf = length(csfIn);

header = cell(Ngdf, 1);
iLowHiCsf = zeros(Ngdf, 1);
g = zeros(Ngdf, 1);
isxcsf = zeros(Ngdf, 1);
isycsf = zeros(Ngdf, 1);
Npatcsf = zeros(Ngdf, 1);
icdef = zeros(Ngdf, 1);

verts = cell(Ngdf, 1);
Ng = cell(Ngdf, 1);
Nk = cell(Ngdf, 1);
knts = cell(Ngdf, 1);

for m = 1:Ngdf
    fid = fopen([sourceFolder '\' csfIn{m} '.csf']);
    
    header{m} = fgetl(fid);
    num = textscan(fid,'%f',1);
    num = num{1};
    iLowHiCsf(m) = num(1);

    fgetl(fid);
    num = textscan(fid,'%f',2);
    num = num{1};
    isxcsf(m) = num(1);
    isycsf(m) = num(2);

    fgetl(fid);
    num = textscan(fid,'%f',3);
    num = num{1};
    Npatcsf(m) = num(1);
    icdef(m) = num(2);
    pszcsf(m) = num(3);
    fgetl(fid);
    
    Ng{m} = zeros(Npatcsf(m), 2);
    Nk{m} = zeros(Npatcsf(m), 2);
    knts{m} = cell(Npatcsf(m), 2);
    verts{m} = cell(Npatcsf(m), 1);
    
    for n = 1:Npatcsf(m)
        num = textscan(fid,'%f',4);
        Nug = num{1}(1);
        Nvg = num{1}(2);
        Kug = num{1}(3);
        Kvg = num{1}(4);
        
        Ng{m}(n,1) = Nug;
        Ng{m}(n,2) = Nvg;
        Nk{m}(n,1) = Kug;
        Nk{m}(n,2) = Kvg;

        Nua = Nug+2*Kug-1;
        Nva = Nvg+2*Kvg-1;
        
        knts{m}(n,1) = textscan(fid,'%f',Nua);
        knts{m}(n,2) = textscan(fid,'%f',Nva);

        Nb = (Nug + Kug - 1)*(Nvg + Kvg - 1);

        num = textscan(fid, '%f', Nb*3);

        verts{m}{n} = reshape(num{1},3,Nb)';
    end
    
    fclose(fid);
end

if ~all(iLowHiCsf == iLowHiCsf(1)) ...
        || ~all(isxcsf == isxcsf(1)) || ~all(isycsf == isycsf(1)) 
    warning('Csf settings are not the same.');
end

for m = 1:Ngdf
    filename = [destFolder '\' csfOut{m} '.csf'];
    fid = fopen(filename, 'wt');

    fprintf(fid, header{m});
    fprintf(fid, '\n%i\tILOWHICSF\n', iLowHiCsf(m));
    fprintf(fid, '%i %i\tISXCSF, ISYCSF\n', isxcsf(m), isycsf(m));
    fprintf(fid, '%i %i %4.2f\tNPATCSF, ICDEF, PSZCSF\n', Npatcsf(m), icdef(m), pszcsf(m)*scale); % The panel size needs adjusting from the default if the body has been scaled!

    for n = 1:Npatcsf(m)
        
        fprintf(fid, '%i %i\n', Ng{m}(n,1), Ng{m}(n,2));
        fprintf(fid, '%i %i\n', Nk{m}(n,1), Nk{m}(n,2));
        
        for ii = 1:2
            for i = 1:length(knts{m}{n,ii})
                fprintf(fid, ' %11.7f', knts{m}{n,ii}(i));
            end
            fprintf(fid, '\n');
        end
        
        Nv = size(verts{m}{n}, 1);
        for ii = 1:Nv
            for i = 1:3
                fprintf(fid, ' %11.4f', scale*(verts{m}{n}(ii,i)) + delta(i)); % The precision here was reduced to sidestep an error where wamit thought the panels were above the free surface
            end
            fprintf(fid, '\n');
        end
        fprintf(fid, '\n');
    end

    fclose(fid);
end

end