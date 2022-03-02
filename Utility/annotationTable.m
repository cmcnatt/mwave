function [] = annotationTable(fig, T, varargin)

TString = evalc('disp(T)');

% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');

% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');

% Output the table using the annotation command.
annotation(fig, 'Textbox', 'String', TString, 'Interpreter', 'Tex', ...
    'FontName', FixedWidth, 'Units', 'Normalized', varargin{:});
end

