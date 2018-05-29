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
classdef IStlGeo < handle
    
    properties (Abstract)
        Vertices;
        Faces;
        Normals;
        Areas;
        Centers;
        DOF;
    end
    
    methods (Abstract)
        SetSetting(stl, val)
        PointsAtSetting(stl)
    end
    
    methods
        function [] = plot(stl, varargin)
            stl.plotFuncs(false, varargin{:});
        end
        
        function [] = surf(stl, varargin)
            stl.plotFuncs(true, varargin{:});
        end
    end
    
    methods (Access = protected)
        function [] = plotFuncs(stl, isSurf, varargin)
            [opts, args] = checkOptions({{'norm'}, {'color', 1}, {'origin'}}, varargin);
            showNorm = opts(1);
            color = [];
            if opts(2)
                color = args{2};
            end
            atOrg = opts(3);
           
            tri = stl.Faces;
            
            if atOrg
                p = stl.Vertices;
                c = stl.Centers;
                n = stl.Normals;
            else
                [p, n, c] = stl.PointsAtSetting;
            end
            
            if ~iscell(tri)
                tri = {tri};
            end
            if ~iscell(p)
                p = {p};
            end
            if ~iscell(n)
                n = {n};
            end
            if ~iscell(c)
                c = {c};
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
            
            for i = 1:length(p)
                if isSurf
                    trisurf(tri{i}, p{i}(:,1), p{i}(:,2), p{i}(:,3), plotOpts{:});
                else
                    trimesh(tri{i}, p{i}(:,1), p{i}(:,2), p{i}(:,3), plotOpts{:});
                end

                if showNorm
                    quiver3(c{i}(:,1), c{i}(:,2), c{i}(:,3), n{i}(:,1), n{i}(:,2), n{i}(:,3))
                end
            end
        end
    end
    
    methods (Static, Access = protected)
                
        function [areas] = triArea(faces, verts)
            % Function to calculate the area of a triangle
            v1 = verts(faces(:,3),:) - verts(faces(:,1),:);
            v2 = verts(faces(:,2),:) - verts(faces(:,1),:);
            av_tmp =  1/2.*(cross(v1,v2));
            areas = sqrt(av_tmp(:,1).^2 + av_tmp(:,2).^2 + av_tmp(:,3).^2);
        end
        
        function [norms] = triNorm(faces, verts)
            TR = triangulation(faces, verts);
            norms = faceNormal(TR);
        end

        function [] = checkNorms(norms)
            % Function to check STL file
            normMag = sqrt(norms(:,1).^2 + norms(:,2).^2 + norms(:,3).^2);
            
            check = (sum(abs(normMag - 1) < 1e-10) + sum(normMag == 0))/length(normMag);
            if check > 1.01 || check < 0.99
                error(['length of normal vectors in stl is not equal to one.'])
            end
        end

        function [cents] = triCenter(faces, verts)
            %Function to caculate the center coordinate of a triangle
            cents = zeros(length(faces),3);
            cents(:,1) = (verts(faces(:,1),1) + verts(faces(:,2),1) + verts(faces(:,3),1))./3;
            cents(:,2) = (verts(faces(:,1),2) + verts(faces(:,2),2) + verts(faces(:,3),2))./3;
            cents(:,3) = (verts(faces(:,1),3) + verts(faces(:,2),3) + verts(faces(:,3),3))./3;
        end
    end
end