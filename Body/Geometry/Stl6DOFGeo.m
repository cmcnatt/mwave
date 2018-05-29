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
classdef Stl6DOFGeo < IStlGeo
    
    properties (Access = private)
        verts;
        faces;
        norms;
        areas;
        cents;
        cg;
        cr;
        pos;
        orient;
    end
    
    properties (Dependent)
        Vertices;
        Faces;
        Normals;
        Areas;
        Centers;
        DOF;
        Cg;
        CenterRot;
        Position;
        Orientation;
        RotationMatrix;
    end
    
    methods
        
        function [stl] = Stl6DOFGeo(fileName)
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
        
        function [] = set.CenterRot(stl, val)
            stl.cr = val;
        end
        function [val] = get.CenterRot(stl)
            val = stl.cr;
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
        
        function [val] = get.DOF(stl)
            val = 6;
        end
        
        function [] = SetSetting(stl, val)
            if length(val) ~= stl.DOF
                error('The value of the setting must be equal to the DOF');
            end
            
            stl.pos = val(1:3);
            stl.orient = val(4:6);
        end
        
        function [val] = get.RotationMatrix(stl)
            % yaw, pitch, roll rotation matrix
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
                val = Rz*Ry*Rx;
            end
        end
        
        function [] = SetStl(stl, verts, faces, cg)
            if nargin < 4
                cg = [0 0 0];
            end
            
            [Nv, ~] = size(verts);
            cgV = ones(Nv, 1)*cg;
            
            stl.verts = verts - cgV;
            stl.faces = faces;
            stl.norms = IStlGeo.triNorm(stl.faces, stl.verts);
            stl.checkNorms(stl.norms);
            stl.areas = IStlGeo.triArea(stl.faces, stl.verts);
            stl.cents = IStlGeo.triCenter(stl.faces, stl.verts);
        end
        
        function [] = Read(stl, fileName, cgFile)
            if nargin < 3
                cgFile = [];
            end
            if isempty(cgFile)
                cgFile = [0 0 0];
            end
            [v, f, n] = importStl(fileName);
            
            [Nv, ~] = size(v);
            cgV = ones(Nv, 1)*cgFile;
            
            stl.verts = v - cgV;
            stl.faces = f;
            stl.norms = n;
            IStlGeo.checkNorms(stl.norms);
            stl.areas = IStlGeo.triArea(stl.faces, stl.verts);
            stl.cents = IStlGeo.triCenter(stl.faces, stl.verts);
        end
        
        function [verts, cents, norms] = PointsAtSetting(stl)
            verts = stl.verts;
            cents = stl.cents;
            norms = stl.norms;
            
            [Nv, ~] = size(verts);
            [Nc, ~] = size(cents);
            
            cg0 = stl.cg;
            if isempty(cg0)
                cg0 = [0 0 0];
            end
            
            cr0 = stl.cr;
            if isempty(cr0);
                cr0 = [0 0 0];
            end
            
            cgV = ones(Nv, 1)*cg0;
            cgC = ones(Nc, 1)*cg0;
            crV = ones(Nv, 1)*cr0;
            crC = ones(Nc, 1)*cr0;
                        
            if ~isempty(stl.orient)
                verts = verts + cgV - crV;
                cents = cents + cgC - crC;
                R = stl.RotationMatrix;
                verts = (R*verts')';
                cents = (R*cents')';
                norms = (R*norms')';
                verts = verts - cgV + crV;
                cents = cents - cgC + crC;
            end
            
            if ~isempty(stl.pos)
                verts = verts + ones(Nv, 1)*stl.pos;
                cents = cents + ones(Nc, 1)*stl.pos;
            end
            
            verts = verts + cgV;
            cents = cents + cgC;
        end
    end
end