function [uR] = rotateAboutPoint(R, point, u)

[N, ~] = size(u);
uR = zeros(size(u));

if iscolumn(point)
    point = point.';
end

for n = 1:N
    uR(n,:) = (R*(u(n,:) - point).').' + point;
end