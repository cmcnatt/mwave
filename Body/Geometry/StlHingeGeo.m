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
classdef StlHingeGeo < IStlGeo
    
    properties (Access = private)
        fwd;
        aft;
        fwdCg;
        aftCg;
        hinLoc;
        hinVec;
        pos;
        orient;
        flex;
    end
    
    properties (Dependent)
        Vertices;
        Faces;
        Normals;
        Areas;
        Centers;
        DOF;
        ForwardBody
        AftBody;
        ForwardCg;
        AftCg;
        HingeLoc;
        HingeVec;
        Position;
        Orientation;
        FlexAngle;
    end
    
    methods
        function [stl] = StlHingeGeo(fwdFileName, aftFileName)
            stl.fwd = Stl6DOFGeo;
            stl.aft = Stl6DOFGeo;
            
            if nargin < 1
                fwdFileName = [];
            end
            
            if nargin < 2
                aftFileName = [];
            end
            
            if ~isempty(fwdFileName)
                stl.fwd.Read(fwdFileName)
            end
            
            if ~isempty(aftFileName)
                stl.aft.Read(aftFileName)
            end
        end
                
        function [] = ReadForward(stl, fileName)
            stl.fwd.Read(fileName);
        end
        
        function [] = ReadAft(stl, fileName)
            stl.aft.Read(fileName);
        end
        
        function [val] = get.Vertices(stl)
            val = {stl.fwd.Vertices, stl.aft.Vertices};
        end
        
        function [val] = get.Faces(stl)
            val = {stl.fwd.Faces, stl.aft.Faces};
        end
        
        function [val] = get.Normals(stl)
            val = {stl.fwd.Normals, stl.aft.Normals};
        end
        
        function [val] = get.Areas(stl)
            val = {stl.fwd.Areas, stl.aft.Areas};
        end
        
        function [val] = get.Centers(stl)
            val = {stl.fwd.Centers, stl.aft.Centers};
        end
        
        function [] = set.ForwardBody(stl, val)
            stl.fwd = val;
        end
        function [val] = get.ForwardBody(stl)
            val = stl.fwd;
        end
        
        function [] = set.AftBody(stl, val)
            stl.aft = val;
        end
        function [val] = get.AftBody(stl)
            val = stl.aft;
        end
        
        function [] = set.ForwardCg(stl, val)
            stl.fwdCg = val;
        end
        function [val] = get.ForwardCg(stl)
            val = stl.fwdCg;
        end
        
        function [] = set.AftCg(stl, val)
            stl.aftCg = val;
        end
        function [val] = get.AftCg(stl)
            val = stl.aftCg;
        end
        
        function [] = set.HingeLoc(stl, val)
            stl.hinLoc = val;
        end
        function [val] = get.HingeLoc(stl)
            val = stl.hinLoc;
        end
        
        function [] = set.HingeVec(stl, val)
            stl.hinVec = val;
        end
        function [val] = get.HingeVec(stl)
            val = stl.hinVec;
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
        
        function [] = set.FlexAngle(stl, val)
            stl.flex = val;
        end
        function [val] = get.FlexAngle(stl)
            val = stl.flex;
        end
        
        function [val] = get.DOF(stl)
            val = 7;
        end
        
        function [] = SetSetting(stl, val)
            if length(val) ~= stl.DOF
                error('The value of the setting must be equal to the DOF');
            end
            
            stl.pos = val(1:3);
            stl.orient = val(4:6);
            stl.flex = val(7);
        end
        
        function [verts, cents, norms] = PointsAtSetting(stl)
            stl.fwd.Cg = stl.fwdCg;
            stl.fwd.CenterRot = stl.fwdCg;
            stl.aft.Cg = stl.aftCg;
            stl.aft.CenterRot = stl.fwdCg;
            stl.fwd.Position = stl.pos;
            stl.aft.Position = stl.pos;
            stl.fwd.Orientation = stl.orient;
            stl.aft.Orientation = stl.orient;
            
            [vf, cf, nf] = stl.fwd.PointsAtSetting;
            [va, ca, na] = stl.aft.PointsAtSetting;
            
            hin0 = stl.hinLoc;
            
            if ~isempty(hin0) && ~isempty(stl.flex)
                                
                rhin = hin0 - stl.fwdCg;
                R = stl.fwd.RotationMatrix;
                rhin = (R*rhin')';
                hin0 = stl.fwdCg + rhin + stl.pos;
                
                if isempty(stl.hinVec)
                    hinVec0 = [0 1 0];
                end
                hinVec0 = (R*hinVec0')';
                
                [Nv, ~] = size(va);
                [Nc, ~] = size(ca);
                
                hinV = ones(Nv, 1)*hin0;
                hinC = ones(Nc, 1)*hin0;
                
                va = va - hinV;
                ca = ca - hinC;
                
                R = stl.hingeRotationMatrix(hinVec0);

                va = (R*va')';
                ca = (R*ca')';
                na = (R*na')';
                
                va = va + hinV;
                ca = ca + hinC;
            end
            
            verts = {vf, va};
            cents = {ca, cf};
            norms = {na, nf};         
        end
    end
    
    methods (Access = private)
        function [R] = hingeRotationMatrix(stl, vect)
            v = vect./sqrt(vect(1)^2 + vect(2)^2 + vect(3)^2);
            
            c = cos(stl.flex);
            s = sin(stl.flex);
            
            x = v(1);
            y = v(2);
            z = v(3);
            
            R = [c + x^2*(1 - c), x*y*(1 - c) - z*s, x*z*(1 - c) + y*s;
                y*x*(1 - c) + z*s, c + y^2*(1 - c), y*z*(1 - c) - x*s;
                z*x*(1 - c) - y*s, z*y*(1 - c) + x*s, c + z^2*(1 - c)];
        end
    end
end