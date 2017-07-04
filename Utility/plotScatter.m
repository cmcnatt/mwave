function [] = plotScatter(x, y, mat, varargin)

[opts, args] = checkOptions({{'xinds', 1}, {'yinds', 1}}, varargin);
     
xinds = 1:length(x);
yinds = 1:length(y);

if opts(1)
    xinds = args{1};
end
if opts(2)
    yinds = args{2};
end

yLab = cell(1,length(y));
for n = 1:length(y)
    yLab{n} = num2str(y(n), '%4.1f');
end

xLab = cell(1,length(x));
for n = 1:length(x)
    xLab{n} = num2str(x(n), '%4.1f');
end

imagesc(mat);
set(gca, 'XAxisLocation','top',...
    'ytick', yinds, 'yticklabel', yLab(yinds),...
    'xtick', xinds, 'xticklabel', xLab(xinds));
fet;

