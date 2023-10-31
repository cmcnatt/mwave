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
classdef StlWave < IStlGeo
    
    properties (Access = private)
        verts;
        faces;
        norms;
        areas;
        cents;
        a;
        omega;
        beta;
        h;
        eta;
        x;
        y;
        time;
        ramp;
    end
    
    properties (Dependent)
        Vertices;
        Faces;
        Normals;
        Areas;
        DOF;
        Centers;
        A;
        Omega;
        Beta;
        H;
        X;
        Y;
        Time;
        Ramp;
    end
    
    methods
        
        function [stl] = StlWave(a, omega, beta, h)
            if nargin > 0
                stl.a = a;
            end
            if nargin > 1
                stl.omega = omega;
            end
            if nargin > 2
                stl.beta = beta;
            end
            if nargin > 3
                stl.h = h;
            end
        end
        
        function [val] = get.Vertices(stl)
            if isempty(stl.verts)
                stl.makeStl;
            end
            val = stl.verts;
        end
        
        function [val] = get.Faces(stl)
            if isempty(stl.faces)
                stl.makeStl;
            end
            val = stl.faces;
        end
        
        function [val] = get.Normals(stl)
            if isempty(stl.norms)
                stl.makeStl;
            end
            val = stl.norms;
        end
        
        function [val] = get.Areas(stl)
            if isempty(stl.areas)
                stl.makeStl;
            end
            val = stl.areas;
        end
        
        function [val] = get.Centers(stl)
            if isempty(stl.cents)
                stl.makeStl;
            end
            val = stl.cents;
        end
        
        function [val] = get.DOF(stl)
            val = 1;
        end
        
        function [] = SetSetting(stl, val)
            if length(val) ~= stl.DOF
                error('The value of the setting must be equal to the DOF');
            end
            
            stl.Time = val;
        end
        
        function [] = set.A(stl, val)
            stl.a = val;
            stl.computeEta;
        end
        function [val] = get.A(stl)
            val = stl.a;
        end
        
        function [] = set.Omega(stl, val)
            stl.omega = val;
            stl.computeEta;
        end
        function [val] = get.Omega(stl)
            val = stl.omega;
        end
        
        function [] = set.Beta(stl, val)
            stl.beta = val;
            stl.computeEta;
        end
        function [val] = get.Beta(stl)
            val = stl.beta;
        end
        
        function [] = set.H(stl, val)
            stl.h = val;
            stl.computeEta;
        end
        function [val] = get.H(stl)
            val = stl.h;
        end
        
        function [] = set.X(stl, val)
            stl.x = val;
            stl.computeEta;
        end
        function [val] = get.X(stl)
            val = stl.x;
        end
        
        function [] = set.Y(stl, val)
            stl.y = val;
            stl.computeEta;
        end
        function [val] = get.Y(stl)
            val = stl.y;
        end
        
        function [] = set.Time(stl, val)
            stl.time = val;
            stl.makeStl;
        end
        function [val] = get.Time(stl)
            val = stl.time;
        end
        
        function [] = set.Ramp(stl, val)
            stl.ramp = val;
        end
        function [val] = get.Ramp(stl)
            val = stl.ramp;
        end
        
        function [verts, cents, norms] = PointsAtSetting(stl)
            
            if isempty(stl.verts) || isempty(stl.cents) || isempty(stl.cents)
                stl.makeStl;
            end
            verts = stl.verts;
            cents = stl.cents;
            norms = stl.norms;
        end
        
        function [maxz, minz] = ComputeMaxMin(stl, time)
            maxz = -inf;
            minz = inf;
            
            for n = 1:length(time)
                stl.Time = time(n);
                verts = stl.PointsAtSetting;
                z = verts(:,3);
                if max(z) > maxz
                    maxz = max(z);
                end
                if min(z) < minz
                    minz = min(z);
                end
            end
        end
    end
        
    methods (Access = private)
        function [] = computeEta(stl)
            if ~isempty(stl.x) && ~isempty(stl.y) && ~isempty(stl.a) ...
                    && ~isempty(stl.omega) && ~isempty(stl.beta) && ~isempty(stl.h)
                
                [X_, Y_] = meshgrid(stl.x, stl.y);
                
                N = length(stl.a);
                if length(stl.omega) ~= N
                    error('N omega must equal N a');
                end
                if length(stl.beta) == 1
                    beta0 = stl.beta*ones(size(stl.a));
                elseif length(stl.beta) ~= N
                    error('N beta must be 1 or equal N a');
                else
                    beta0 = stl.beta;
                end
                
                k = solveForK(stl.omega, stl.h);
                
                stl.eta = zeros(length(stl.y), length(stl.x), N);
                for n = 1:N
                    stl.eta(:,:,n) = stl.a(n)*exp(-1i*(k(n)*cos(beta0(n))*X_ + k(n)*sin(beta0(n))*Y_));
                end   
            end
        end
        
        function [] = makeStl(stl)
            
            if isempty(stl.eta)
                stl.computeEta;
            end
            if ~isempty(stl.x) && ~isempty(stl.y) && ~isempty(stl.a) ...
                    && ~isempty(stl.omega) && ~isempty(stl.beta) && ~isempty(stl.h)
                
                if isempty(stl.time)
                    stl.time = 0;
                end
                
                N = length(stl.a);
                etat = zeros(length(stl.x), length(stl.y));
                if isempty(stl.ramp)
                    ram = 1;
                else
                    ram = stl.ramp;
                end
                for n = 1:N
                    etan = reshape(stl.eta(:,:,n), [length(stl.x) length(stl.y)]);
                    etat = etat + ram*real(etan*exp(1i*stl.omega(n)*stl.time));
                end

                [X_, Y_] = meshgrid(stl.x, stl.y);

                N = numel(X_);
                x_ = reshape(X_, [N, 1]);
                y_ = reshape(Y_, [N, 1]);
                z_ = reshape(etat, [N, 1]);

                stl.faces = delaunay(x_, y_);
                stl.verts = [x_, y_, z_];
                stl.norms = IStlGeo.triNorm(stl.faces, stl.verts);
                stl.areas = IStlGeo.triArea(stl.faces, stl.verts);
                stl.cents = IStlGeo.triCenter(stl.faces, stl.verts);
            end
            
        end
    end
end