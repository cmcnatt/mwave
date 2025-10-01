function Wamit_writeCsfHi(folderPath, geoName, panelGeoCsf, pszcsf)

filename = fullfile(folderPath, [geoName, '.csf']);
fileID = fopen(filename, 'wt');

ilowhicsf = 1; % This function is for higher order CSF file creation.
icdef = 0; % 0 for flat quadrilateral panels

fprintf(fileID, ['Model ' geoName ', created: ' date '\n']);
fprintf(fileID, '%i \n', ilowhicsf);
fprintf(fileID, '%i %i \n', panelGeoCsf.Xsymmetry, panelGeoCsf.Ysymmetry);
fprintf(fileID, '%i %i %8.4f\n', panelGeoCsf.Count, icdef, pszcsf);

pans = panelGeoCsf.Panels;

for n = 1:panelGeoCsf.Count
    pan = pans(n);
    
    verts = pan.Vertices;
    for m = 1:4
        fprintf(fileID, '\t%8.4f\t%8.4f\t%8.4f\n', verts(m,1), verts(m,2), verts(m,3));
    end
end

fclose(fileID);

end

