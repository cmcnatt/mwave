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
function [nondim] = nondimMats(types, iheave, M, C, varargin)

[opts, args] = checkOptions({{'A', 1}, {'B', 1}, {'D', 1}, {'K', 1},...
    {'f', 1}, {'omega', 1}, {'lam', 1}, {'L', 1}}, varargin);

dof = length(types);

checkSize(dof, M);
checkSize(dof, C);

nondim.M = zeros(dof, dof);
nondim.C = zeros(dof, dof);

nf = 1;

if (opts(1))
    noA = false;
    A = args{1};
    checkSize(dof, A);
    if (2 == ndims(A))
        nondim.A = zeros(dof, dof);
    else
        [nf, ~, ~] = size(A);
        nondim.A = zeros(nf, dof, dof);
    end
else
    noA = true;
end

if (opts(2))
    noB = false;
    B = args{2};
    checkSize(dof, B);
    if (2 == ndims(B))
        nondim.B = zeros(dof, dof);
    else
        [nf, ~, ~] = size(B);
        nondim.B = zeros(nf, dof, dof);
    end
else
    noB = true;
end

if (opts(3))
    noD = false;
    D = args{3};
    checkSize(dof, D);
    nondim.D = zeros(dof, dof);
else
    noD = true;
end

if (opts(4))
    noK = false;
    K = args{4};
    checkSize(dof, K);
    nondim.K = zeros(dof, dof);
else
    noK = true;
end

if (opts(5))
    nof = false;
    f = args{5};
    if (dof ~= length(f))
        error('The length of the force vector must be equal to the DOF');
    end
    nondim.f = zeros(dof, 1);
else
    nof = true;
end

if (opts(6))
    noOmega = false;
    omega = args{6};
else
    noOmega = true;
end

if (opts(7))
    noLam = false;
    lam = args{7};
else
    noLam = true;
end

if (opts(8))
    noL = false;
    L = args{8};
else
    noL = true;
end

% Nondimensionalizing variables
m = M(iheave, iheave);
c = C(iheave, iheave);
g = IWaves.G;

for n = 1:dof
    for o = 1:dof
        if ((1 == types(n)) && (1 == types(o)))
            % trans-trans
            
            nondim.M(n, o) = M(n, o)/m;
            if (~noA)
                if (1 == nf)
                    nondim.A(n, o) = A(n, o)/m;
                else
                    nondim.A(:, n, o) = A(:, n, o)./m;
                end
            end
            
            nondim.C(n, o) = C(n, o)/c;
            if (~noK)
                nondim.K(n, o) = K(n, o)/c;
            end
        elseif (((1 == types(n)) && (2 == types(o)))...
                || ((2 == types(n)) && (1 == types(o))))
            % trans-rot
            
            nondim.M(n, o) = c/(m^2*g)*M(n, o);
            if (~noA)
                if (1 == nf)
                    nondim.A(n, o) = c/(m^2*g)*A(n, o);
                else
                    nondim.A(:, n, o) = c/(m^2*g)*A(:, n, o);
                end
            end
            
            nondim.C(n, o) = (1/(m*g))*C(n, o);
            if (~noK)
                nondim.K(n, o) = (1/(m*g))*K(n, o);
            end
        elseif ((2 == types(n)) && (2 == types(o)))
            % rot-rot
            
            nondim.M(n, o) = c^2/(m^3*g^2)*M(n, o);
            if (~noA)
                if (1 == nf)
                    nondim.A(n, o) = c^2/(m^3*g^2)*A(n, o);
                else
                    nondim.A(:, n, o) = c^2/(m^3*g^2)*A(:, n, o);
                end
            end
            
            nondim.C(n, o) = (c/(m^2*g^2))*C(n, o);
            if (~noK)
                nondim.K(n, o) = (c/(m^2*g^2))*K(n, o);
            end
        else
            error(['types should be ''1'' to indicate translational ',...
                'modes or ''2'' to indicate rotational modes']);
        end
    end
end

if (~noOmega)
    nondim.omega = sqrt(m/c)*omega;
end

if (~noLam)
    nondim.lam = c/(m*g)*lam;
end

if (~noL)
    nondim.L = c/(m*g)*L;
end

end

function [] = checkSize(dof, mat)
    dim2 = dof;
    switch ndims(mat)
        case 1
            M = length(mat);
            N = 1;
            dim2 = 1;
        case 2
            [M, N] = size(mat);
        case 3
            [~, M, N] = size(mat);
    end
    
    if ((M ~= dof) || (N ~= dim2))
        error('The size of each matrix must be equal to the DOF');
    end
end