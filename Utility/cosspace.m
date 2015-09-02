function [x] = cosspace(x1, x2, varargin)

if (x2 == x1)
    error('cosspace: The bounds (x1, x2) cannot be the same value.');
elseif (x2 < x1)
    neg = true;
    x1p = x1;
    x1 = x2;
    x2 = x1p;
else
    neg = false;
end

if (~isempty(varargin))
    N = varargin{1};
else
    N = 101;
end

dtheta = linspace(-pi,0,N);

xlen = x2 - x1;

x = xlen/2*(ones(size(dtheta)) + cos(dtheta)) + x1*ones(size(dtheta));

if (neg)
    x = flip(x);
end