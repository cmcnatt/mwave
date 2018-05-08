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
classdef StlGeo < handle
    
    properties (Access = private)
        verts;
        faces;
        norms;
        areas;
        cents;
        cg;
        pos;
        orient;
    end
    
    properties (Dependent)
        Vertices;
        Faces;
        Normals;
        Areas;
        Centers;
        Cg;
        Position;
        Orientation;
        RotationMatrix;
    end
    
    methods
        
        function [stl] = StlGeo(fileName)
            if nargin < 1
                fileName = [];
            end
            
            if ~isempty(fileName)
                stl.Read(fileName)
            end
        end
        
        function [val] = get.Vertices(stl)
            val = stl.verts;
        end
        
        function [val] = get.Faces(stl)
            val = stl.faces;
        end
        
        function [val] = get.Normals(stl)
            val = stl.norms;
        end
        
        function [val] = get.Areas(stl)
            val = stl.areas;
        end
        
        function [val] = get.Centers(stl)
            val = stl.cents;
        end
        
        function [] = set.Cg(stl, val)
            stl.cg = val;
        end
        function [val] = get.Cg(stl)
            val = stl.cg;
        end
        
        function [] = set.Position(stl, val)
            stl.pos = val;
        end
        function [val] = get.Position(stl)
            val = stl.pos;
        end
        
        function [] = set.Orientation(stl, val)
            stl.orient = val;
        end
        function [val] = get.Orientation(stl)
            val = stl.orient;
        end
        
        function [val] = get.RotationMatrix(stl)
            a = stl.orient;
            if isempty(a)
                val = [];
            else
                Rx = [1 0 0;
                    0 cos(a(1)) -sin(a(1));
                    0 sin(a(1)) cos(a(1))];
                Ry = [cos(a(2)) 0 sin(a(2));
                    0 1 0;
                    -sin(a(2)) 0 cos(a(2))];
                Rz = [cos(a(3)) -sin(a(3)) 0;
                    sin(a(3)) cos(a(3)) 0;
                    0 0 1];
                val = Rx*Ry*Rz;
            end
        end
        
        function [] = Read(stl, fileName) 
            [stl.verts, stl.faces, stl.norms] = importStl(fileName);
            stl.checkStl();
            stl.triArea();
            stl.triCenter();
        end
        
        function [verts, cents, norms] = PointsAtPosOrient(stl)
            verts = stl.verts;
            cents = stl.cents;
            norms = stl.norms;
            
            [Nv, ~] = size(verts);
            [Nc, ~] = size(cents);
            
            cg0 = stl.cg;
            if isempty(cg0)
                cg0 = [0 0 0];
            end
            
            cgV = ones(Nv, 1)*cg0;
            cgC = ones(Nc, 1)*cg0;
            
            verts = verts - cgV;
            cents = cents - cgC;
            
            if ~isempty(stl.orient)
                R = stl.RotationMatrix;
                verts = (R*verts')';
                cents = (R*cents')';
                norms = (R*norms')';
            end
            
            if ~isempty(stl.pos)
                verts = verts + ones(Nv, 1)*stl.pos;
                cents = cents + ones(Nc, 1)*stl.pos;
            end
            
            verts = verts + cgV;
            cents = cents + cgC;
        end
                
        function [] = plot(stl, varargin)
            stl.plotFuncs(false, varargin{:});
        end
        
        function [] = surf(stl, varargin)
            stl.plotFuncs(true, varargin{:});
        end
    end
    
    methods (Access = private)
        function [] = plotFuncs(stl, isSurf, varargin)
            [opts, args] = checkOptions({{'norm'}, {'color', 1}, {'origin'}}, varargin);
            showNorm = opts(1);
            color = [];
            if opts(2)
                color = args{2};
            end
            atOrg = opts(3);
           
            tri = stl.faces;
            
            if atOrg
                p = stl.verts;
                c = stl.cents;
                n = stl.norms;
            else
                [p, n, c] = stl.PointsAtPosOrient;
            end
            
            hold on 
            
            plotOpts = {};
            if ~isempty(color)
                if isSurf
                    plotOpts = {'FaceColor', color};
                else
                    plotOpts = {'EdgeColor', color};
                end
            end
            
            if isSurf
                trisurf(tri, p(:,1), p(:,2), p(:,3), plotOpts{:});
            else
                trimesh(tri, p(:,1), p(:,2), p(:,3), plotOpts{:});
            end
            
            if showNorm
                quiver3(c(:,1), c(:,2), c(:,3), n(:,1), n(:,2), n(:,3))
            end
        end
                
        function [] = triArea(stl)
            % Function to calculate the area of a triangle
            v1 = stl.verts(stl.faces(:,3),:)-stl.verts(stl.faces(:,1),:);
            v2 = stl.verts(stl.faces(:,2),:)-stl.verts(stl.faces(:,1),:);
            av_tmp =  1/2.*(cross(v1,v2));
            area_mag = sqrt(av_tmp(:,1).^2 + av_tmp(:,2).^2 + av_tmp(:,3).^2);
            stl.areas = area_mag;
        end

        function [] = checkStl(stl)
            % Function to check STL file
            tnorm = stl.norms;
            norm_mag = sqrt(tnorm(:,1).^2 + tnorm(:,2).^2 + tnorm(:,3).^2);
            check = sum(norm_mag)/length(norm_mag);
            if check > 1.01 || check < 0.99
                error(['length of normal vectors in stl is not equal to one.'])
            end
        end

        function [] = triCenter(stl)
            %Function to caculate the center coordinate of a triangle
            points = stl.verts;
            c = zeros(length(stl.faces),3);
            c(:,1) = (points(stl.faces(:,1),1)+points(stl.faces(:,2),1)+points(stl.faces(:,3),1))./3;
            c(:,2) = (points(stl.faces(:,1),2)+points(stl.faces(:,2),2)+points(stl.faces(:,3),2))./3;
            c(:,3) = (points(stl.faces(:,1),3)+points(stl.faces(:,2),3)+points(stl.faces(:,3),3))./3;
            stl.cents = c;
        end
    end
end