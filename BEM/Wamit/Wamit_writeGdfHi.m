function Wamit_writeGdfHi(folderPath, geoName, panelGeo, excludeInt, onlyWets)

if nargin < 4
    excludeInt = false;
end

if nargin < 5
    onlyWets = false;
end

filename = fullfile(folderPath, [geoName, '.gdf']);
fileID = fopen(filename, 'wt');

ulen = 1;
g = 9.806650;
igdef = 0;

if onlyWets
    N = panelGeo.Count;
else
    N = sum(panelGeo.IsWets);
end

Nint = sum(panelGeo.IsInteriors);
if excludeInt
    count = N - Nint;
else
    count = N;
end

fprintf(fileID, ['Model ' geoName ', created: ' date '\n']);
fprintf(fileID, '%8.4f %8.4f\n', ulen, g);
fprintf(fileID, '%i %i \n', panelGeo.Xsymmetry, panelGeo.Ysymmetry);
fprintf(fileID, '%i %i\n', count, igdef);

pans = panelGeo.Panels;

for n = 1:panelGeo.Count
    pan = pans(n);

    panOk = true;
    if (onlyWets && ~pan.IsWet)
        panOk = false;
    end
    if (pan.IsInterior && excludeInt)
        panOk = false;
    end

    if panOk
        verts = pan.Vertices;
        for m = 1:4
            fprintf(fileID, '\t%8.4f\t%8.4f\t%8.4f\n', verts(m,1), verts(m,2), verts(m,3));
        end
    end
end

fclose(fileID);

end

