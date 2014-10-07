%{ 
mwave - A water wave and wave energy converter computation package 
Copyright (C) 2014  Cameron McNatt

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Contributors:
    C. McNatt
%}
function [As, Ar, N0, H0pk0a] = cylinderWaveCoeffs(radius, draft, depth, T, M, N)
% function is an attempt to analytically compute the circular scattering and radiated
% wave coefficients for a cylinder.
%
% Has not been successfully completed.

a = radius;
d = depth;
h = depth - draft;
omega = 2*pi./T;

k = WaveComp.SolveForK(omega, depth, WaveComp.G, 'Evanescent', M);
k0 = k(1);
kq = k(2:end);

N0 = 1/2*(1 + sinh(2*k0*d)/(2*k0*d));
Nq = 1/2*(1 + sin(2*kq*d)./(2*kq*d));

% Scattering Coeffs
As = zeros(M + 1, 2*N + 1);

% n = 0;
U0s = makeUns(0, a, h, d, k);

[E0sq, H0pk0a, K0pkqa] = makeEnsq(0, a, h, k, N0, Nq);

G0qs = makeGnqs(0, a, h, d, k, N0, Nq);

LHS = eye(M + 1) + G0qs*E0sq;
RHS = G0qs*U0s;
D0q = LHS\RHS;

As(1, N + 1) = 1i*cosh(k0*d)/(sqrt(N0)*H0pk0a)*D0q(1);
As(2:end, N + 1) = 1i./(sqrt(Nq.').*K0pkqa.').*D0q(2:end);

for n = 1:N
    %-n
    Unns = makeUns(-n, a, h, d, k);

    [Ennsq Hnnpk0a Knnpkqa] = makeEnsq(-n, a, h, k, N0, Nq);

    Gnnqs = makeGnqs(-n, a, h, d, k, N0, Nq);

    LHS = eye(M + 1) + Gnnqs*Ennsq;
    RHS = Gnnqs*Unns;
    Dnnq = LHS\RHS;

    As(1, N + 1 - n) = 1i*cosh(k0*d)/(sqrt(N0)*Hnnpk0a)*Dnnq(1);
    As(2:end, N + 1 - n) = 1i./(sqrt(Nq.').*Knnpkqa.').*Dnnq(2:end);
    
    %n
    Uns = makeUns(n, a, h, d, k);

    [Ensq Hnpk0a Knpkqa] = makeEnsq(n, a, h, k, N0, Nq);

    Gnqs = makeGnqs(n, a, h, d, k, N0, Nq);

    LHS = eye(M + 1) + Gnqs*Ensq;
    RHS = Gnqs*Uns;
    Dnq = LHS\RHS;

    As(1, N + 1 + n) = 1i*cosh(k0*d)/(sqrt(N0)*Hnpk0a)*Dnq(1);
    As(2:end, N + 1 + n) = 1i./(sqrt(Nq.').*Knpkqa.').*Dnq(2:end);
end


% Radiation Coeffs
Qs = makeQs(a, h, omega, M);

Sq = makeSq(a, h, d, omega, k, N0, Nq);

LHS = eye(M + 1) + G0qs*E0sq;
RHS = G0qs*Qs + Sq;
Dq = LHS\RHS;

Ar = zeros(M + 1, 1);
Ar(1) = 1i*cosh(k0*d)/(sqrt(N0)*H0pk0a)*Dq(1);
Ar(2:end) = 1i./(sqrt(Nq.').*K0pkqa.').*Dq(2:end);

end

function [Ensq Hnpk0a Knpkqa] = makeEnsq(n, a, h, k, N0, Nq)

k0 = k(1);
kq = k(2:end);
M = length(kq);
s = 0:M;

Hnk0a = besselh(n, 2, k0*a);
Hnpk0a = 1/2*(besselh(n - 1, 2, k0*a) - besselh(n + 1, 2, k0*a));

Knkqa = besselk(n, kq*a);
Knpkqa = 1/2*(besselk(n - 1, kq*a) + besselk(n + 1, kq*a));

Ens0 = -2/h*Hnk0a/Hnpk0a*(h^2*k0*(-1).^s.*sinh(k0*h))./(sqrt(N0)*(s.^2*pi^2 + k0^2*h^2));

Kng = meshgrid(Knkqa./Knpkqa, ones(1, M + 1));
kqg = meshgrid(kq, ones(1, M + 1));
sg = meshgrid(s, ones(1, M))';
Nqg = meshgrid(Nq, ones(1, M + 1));

Ensq1 = -2/h*Kng.*(h^2*kqg.*(-1).^sg.*sin(h*kqg))./(sqrt(Nqg).*(-sg.^2*pi^2 + kqg.^2*h^2));

Ensq = zeros(M + 1, M + 1);
Ensq(:,1) = Ens0.';
Ensq(:,2:end) = Ensq1;

end

function [Gnqs] = makeGnqs(n, a, h, d, k, N0, Nq)

k0 = k(1);
kq = k(2:end);
M = length(kq);
s = 1:M;

Gn00 = abs(n)*sinh(k0*h)/(2*a*k0^2*d*sqrt(N0));
Gnq0 = abs(n)*sin(kq*h)./(2*a*kq.^2*d.*sqrt(Nq));

Ins = besseli(n, pi*a/h*s);
Inps = 1/2*(besseli(n - 1, pi*a/h*s) + besseli(n + 1, pi*a/h*s));

Gn0s = Inps./Ins.*(s*pi*h.*(-1).^s*sinh(k0*h))./((s.^2*pi^2 + k0^2*h^2)*d*sqrt(N0));

Ing = meshgrid(Inps./Ins, zeros(1, M));
kqg = meshgrid(kq, ones(1, M))';
sg = meshgrid(s, ones(1, M));
Nqg = meshgrid(Nq, ones(1, M))';

Gnqs1 = Ing.*(sg*pi*h.*(-1).^sg.*sin(kqg*h))./((-sg.^2*pi^2 + kqg.^2*h^2)*d.*sqrt(Nqg));

Gnqs = zeros(M + 1, M + 1);
Gnqs(1,1) = Gn00;
Gnqs(2:end,1) = Gnq0.';
Gnqs(1,2:end) = Gn0s;
Gnqs(2:end,2:end) = Gnqs1;

end

function [Uns] = makeUns(n, a, h, d, k)

M = length(k) - 1;
k0 = k(1);
s = 0:M;

Hnpk0a = 1/2*(besselh(n - 1, 2, k0*a) - besselh(n + 1, 2, k0*a));

Uns = 4*1i^(n + 1)*(-1).^s*h*sinh(k0*h)./(pi*a*(s.^2*pi^2 + k0^2*h^2)*Hnpk0a*cosh(k0*d));
Uns = Uns.';

end

function [Qs] = makeQs(a, h, omega, M)

s = 1:M;
g = WaveComp.G;

Q0 = 1i*omega^2/(g*h)*(h^2/3 - a^2/2);
Qs1 = 2*1i*omega^2*h*(-1).^s./(g*s.^2*pi^2);

Qs = zeros(M + 1, 1);
Qs(1) = Q0;

Qs(2:end) = Qs1.';

end

function [Sq] = makeSq(a, h, d, omega, k, N0, Nq)

g = WaveComp.G;
k0 = k(1);
kq = k(2:end);
M = length(kq);

S0 = 1i*a*omega^2/(2*k0^2*d*g*h*sqrt(N0))*sinh(k0*h);
Sq1 = 1i*a*omega^2./(2*kq.^2*d*g*h.*sqrt(Nq)).*sin(kq*h);

Sq = zeros(M + 1, 1);
Sq(1) = S0;
Sq(2:end) = Sq1.';

end

