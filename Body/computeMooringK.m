function [K] = computeMooringK(cg, pos, ks, anghs, angvs)
% computes spring damping matrix for a 6DOF body

nlin = max([length(ks), length(anghs), length(angvs)]);

if length(ks) == 1
    ks = ks*ones(1, nlin);
elseif length(ks) ~= nlin
    error('ks must be a scalar or a vector of length of the number of lines');
end

if length(anghs) == 1
    anghs = anghs*ones(1, nlin);
elseif length(anghs) ~= nlin
    error('anghs must be a scalar or a vector of length of the number of lines');
end

if length(angvs) == 1
    angvs = angvs*ones(1, nlin);
elseif length(angvs) ~= nlin
    error('angvs must be a scalar or a vector of length of the number of lines');
end

npnt = size(pos,1);

if (npnt ~= 1)
    if (npnt ~= nlin)
        error('Must be a single mooring point on the body or a number equal to the number of mooring lines.');
    end
    pos0 = pos;
else
    pos0 = zeros(nlin,3);
    for n = 1:nlin
        pos0(n,:) = pos;
    end
end

K = zeros(6,6);
k0 = zeros(1,3);

for n = 1:nlin
    k0(1) = ks(n)*cos(angvs(n))*abs(cos(anghs(n)));
    k0(2) = ks(n)*cos(angvs(n))*abs(sin(anghs(n)));
    k0(3) = ks(n)*sin(angvs(n));
    
    r = pos0(n,:) - cg;
    rx = skewMat(r);
    Il = zeros(6,3);
    Il(1:3,1:3) = eye(3);
    Il(4:6,1:3) = rx;

    Ir = zeros(3,6);
    Ir(1:3,1:3) = eye(3);
    Ir(1:3,4:6) = -rx;

    for m = 1:3
        Ir(m,:) = k0(m)*Ir(m,:);
    end

    Kn = Il*Ir;
    K = K + Kn;
end

%K = round(10^9*K)/10^9;



