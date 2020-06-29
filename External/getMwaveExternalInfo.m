function info = getMwaveExternalInfo

fileName = [mwavePath '\External\mwave_External_Info.xlsx'];

opts = spreadsheetImportOptions('NumVariables', 5);
opts.VariableNamesRange = 'A1';
opts.DataRange = 'A2';

info = readtable(fileName, opts);

end