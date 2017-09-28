function [] = drawCoor(ax, loc, R, varargin)

[opts, args] = checkOptions({{'arrowLength', 1}, {'arrowWidth', 1}, ...
    {'Color', 1}, {'linestyle', 1}, {'headsize', 1}}, varargin);

arrLen = 1;
if opts(1)
    arrLen = args{1};
end

arrWid = 0.8;
if opts(2)
    arrWid = args{2};
end

cols = zeros(3,3);
if opts(3)
    colIn = args{3};
    [N, ~] = size(colIn);
    if N == 1
        cols = ones(3,1)*colIn;
    else
        for n = 1:N
            cols(n,:) = colIn(n,:);
        end
    end
else
    for n = 1:3
        cols(n,:) = MColor.Colors(n);
    end
end

if opts(4)
    linStyIn = args{4};
    if ~iscell(linStyIn)
        linSty{1} = linStyIn;
        linSty{2} = linStyIn;
        linSty{3} = linStyIn;
    else
        linSty = linStyIn;
    end
else
%     linSty{1} = '-';
%     linSty{2} = '--';
%     linSty{3} = ':';
    linSty{1} = '-';
    linSty{2} = '-';
    linSty{3} = '-';
end

headsize = 0.5;
if opts(5)
    headsize = args{5};
end

is2d = (length(loc) == 2);

if is2d
    for n = 1:2
        quiver(ax, loc(1), loc(2), R(1,n), R(2,n), arrLen, ...
            'Color', cols(n,:), 'linewidth', arrWid, 'linestyle', linSty{n}, ...
            'MaxHeadSize', headsize);
        hold on;
    end
else
    for n = 1:3
        quiver3(ax, loc(1), loc(2), loc(3), R(1,n), R(2,n), R(3,n), arrLen,...
            'Color', cols(n,:), 'linewidth', arrWid, 'linestyle', linSty{n}, ...
            'MaxHeadSize', headsize);
        hold on;
    end
end
    
