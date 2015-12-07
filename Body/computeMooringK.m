function [K] = computeMooringK(cg, pos, ks, anghs, angvs, varargin)

opts = checkOptions({'NoYaw'}, varargin);
noyaw = opts(1);

nlin = length(ks);

npnt = size(pos,1);

if (npnt ~= 1)
    if (npnt ~= nlin)
        error('Must be a single mooring point on the body or a number equal to the number of mooring lines.');
    end
    poss = pos;
else
    poss = zeros(nlin,3);
    for n = 1:nlin
        poss(n,:) = pos;
    end
end


ksu = 0;
ksw = 0;

khe = 0;
kr = 0;
kp = 0;
ky = 0;

for n = 1:nlin
    
    ksu1 = ks(n)*cos(angvs(n))*abs(cos(anghs(n)));
    ksw1 = ks(n)*cos(angvs(n))*abs(sin(anghs(n)));
    khe1 = ks(n)*sin(angvs(n));
    
    ksu = ksu + ksu1;
    ksw = ksw + ksw1;
    khe = khe + khe1;
    
    r = abs(poss(n,:) - cg);
    
    kr = kr + r(2)*khe1 + r(3)*ksw1;
    kp = kp + r(1)*khe1 + r(3)*ksu1;
    ky = ky + r(1)*ksw1 + r(2)*ksu1;
end

K = zeros(6,6);
K(1,1) = ksu;
K(2,2) = ksw;
K(3,3) = khe;
K(4,4) = kr;
K(5,5) = kp;
if (~noyaw)
    K(6,6) = ky;
end

K = round(10^9*K)/10^9;



