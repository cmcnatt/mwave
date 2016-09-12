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

kph = 0;
kps = 0;

krr = 0;
kpp = 0;
kyy = 0;

for n = 1:nlin
    
    ksu_n = ks(n)*cos(angvs(n))*abs(cos(anghs(n)));
    ksw_n = ks(n)*cos(angvs(n))*abs(sin(anghs(n)));
    khe_n = ks(n)*sin(angvs(n));
    
    ksu = ksu + ksu_n;
    ksw = ksw + ksw_n;
    khe = khe + khe_n;
    
    r = poss(n,:) - cg;
    
    kps = kps - r(3)*ksu_n;
    kph = kph + r(1)*khe_n;
    
    krr = krr + r(3)^2*ksw_n + r(2)^2*khe_n;
    kpp = kpp + r(3)^2*ksu_n + r(1)^2*khe_n;
    kyy = kyy + r(2)^2*ksu_n + r(1)^2*ksw_n;
end

K = zeros(6,6);
K(1,1) = ksu;
K(2,2) = ksw;
K(3,3) = khe;
K(1,5) = kps;
K(5,1) = kps;
K(3,5) = kph;
K(5,3) = kph;
K(4,4) = krr;
K(5,5) = kpp;
if (~noyaw)
    K(6,6) = kyy;
end

K = round(10^9*K)/10^9;



