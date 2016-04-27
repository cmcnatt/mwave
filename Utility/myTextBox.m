function [] = myTextBox(pos, varargin)

nargs = length(varargin);
for n = 1:nargs
    text(n) = varargin(n);
end

annotation('textbox', [pos(1), pos(2), 0.1 0.1],...
            'String', text,...
            'LineStyle', 'none');