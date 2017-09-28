function [uT] = transformPoints(R, r0, u)

if isempty(R)
    R = eye(3);
end

if isempty(r0)
    r0 = [0 0 0];
end

[N, ~] = size(u);
uT = zeros(size(u));

for n = 1:N
    uT(n,:) = (R*u(n,:).').' + r0;
end