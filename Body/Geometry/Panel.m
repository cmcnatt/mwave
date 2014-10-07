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
classdef Panel < handle
    
    properties (Access = private)
        vertices;
        normal;
        area;
        centroid;
        value;
    end
    
    properties (Dependent)
        Vertices;
        Normal;
        Area;
        Centroid;
        Value;
    end
    
    methods
        % Constructor
        function [pan] = Panel(verts)
            if (nargin == 1)
                pass = Panel.checkVertices(verts);
                if (pass)
                    pan.vertices = verts;
                end
                pan.normal = [];
                pan.area = [];
                pan.centroid = [];
                pan.value = [];
            end
        end
        
        function [verts] = get.Vertices(pan)
            verts = pan.vertices;
        end
                
        function [norm] = get.Normal(pan)
            if (isempty(pan.normal))
                pan.computeQuantities;
            end
            
            norm = pan.normal;
        end
        
        function [ar] = get.Area(pan)
            if (isempty(pan.area))
                pan.computeQuantities;
            end
            
            ar = pan.area;
        end
        
        function [cent] = get.Centroid(pan)
            if (isempty(pan.centroid))
                pan.computeQuantities;
            end
            
            cent = pan.centroid;
        end
        
        function [val] = get.Value(pan)
            val = pan.value;
        end
        function [pan] = set.Value(pan, val)
            if (~isempty(val))
                if (~isscalar(val))
                    error('Panel value must be a scalar');
                end
            end
            pan.value = val;
        end
        
        function [] = Translate(pan, vector)
            if (length(vector) ~= 3)
                error('Transtaltional vector must be Nx1 or Nx3');
            end
            n = size(vector);
            if (n == 3)
                vector = vector';
            end
            for n = 1:4
                pan.vertices(n,:) = pan.vertices(n,:) + vector;
            end
            
            pan.centroid = [];
        end
        
        function [] = Rotate(pan, varargin)
            
            [opts, args] = checkOptions({{'RotMat', 1}, {'AxisAngle', 2}, {'Origin', 1}}, varargin);
            
            if (~xor(opts(1), opts(2)))
                error('Input must be either a rotation matrix or an axis and angle');
            end
            
            if (opts(1))
                R = args{1};
            else
                ax = args{2}{1};
                angle = args{2}{2};
                R = createRotationMatrix(ax, angle);
            end
            
            if (opts(3))
                transOrg = true;
                vector = args{3};
                if(iscolumn(vector))
                    vector = vector.';
                end
            else
                transOrg = false;
            end
            
            for n = 1:4
                if (transOrg)
                    v = pan.vertices(n,:) - vector;
                    Rv = R*v.';
                    pan.vertices(n,:) = Rv.' + vector;
                else
                    v = R*pan.vertices(n,:).';
                    pan.vertices(n,:) = v.';
                end
            end
            
            pan.centroid = [];
        end
                
        function [] = plot(pan, varargin)
            xsym = false;
            ysym = false;
            showN = false;
            for n = 1:length(varargin)
                if (strcmp(varargin{n}, 'xsym'))
                    xsym = true;
                end
                if (strcmp(varargin{n}, 'ysym'))
                    ysym = true;
                end
                if (strcmp(varargin{n}, 'ShowNorm'))
                    showN = true;
                end
            end
            pnts = zeros(5,3);
            for n = 1:3
                pnts(:,n) = [pan.vertices(:,n); pan.vertices(1,n)];
            end
            plot3(pnts(:,1), pnts(:,2), pnts(:,3));
            
            if (showN)
                hold on;
                cent = pan.Centroid;
                norm = pan.Normal;
                area = pan.Area;
                quiver3(cent(1), cent(2), cent(3), norm(1), norm(2), norm(3), 10*area);
            end
            
            if (xsym && ~ysym)
                hold on;
                plot3(-pnts(:,1), pnts(:,2), pnts(:,3));
            elseif (~xsym && ysym)
                hold on;
                plot3(pnts(:,1), -pnts(:,2), pnts(:,3));
            elseif (xsym && ysym)
                hold on;
                plot3(pnts(:,1), -pnts(:,2), pnts(:,3));
                plot3(-pnts(:,1), pnts(:,2), pnts(:,3));
                plot3(-pnts(:,1), -pnts(:,2), pnts(:,3));
            end
        end
        
        function [] = surf(pan, varargin)
            xsym = false;
            ysym = false;
            showN = false;
            for n = 1:length(varargin)
                if (strcmp(varargin{n}, 'xsym'))
                    xsym = true;
                end
                if (strcmp(varargin{n}, 'ysym'))
                    ysym = true;
                end
                if (strcmp(varargin{n}, 'ShowNorm'))
                    showN = true;
                end
            end
            
            x = pan.vertices(:,1);
            y = pan.vertices(:,2);
            z = pan.vertices(:,3);
            
            X = zeros(2,2);
            Y = zeros(2,2);
            Z = zeros(2,2);
            
            X(1,1) = x(1);
            X(2,1) = x(2);
            X(2,2) = x(3);
            X(1,2) = x(4);
            
            Y(1,1) = y(1);
            Y(2,1) = y(2);
            Y(2,2) = y(3);
            Y(1,2) = y(4);
            
            Z(1,1) = z(1);
            Z(2,1) = z(2);
            Z(2,2) = z(3);
            Z(1,2) = z(4);
            
            if (~isempty(pan.value))
                C = pan.value*ones(2,2);
                surf(X, Y, Z, C);
            else
                surf(X, Y, Z);
            end
            
            if (showN)
                hold on;
                cent = pan.Centroid;
                norm = pan.Normal;
                area = pan.Area;
                quiver3(cent(1), cent(2), cent(3), norm(1), norm(2), norm(3), 10*area);
            end
            
            if (xsym && ~ysym)
                hold on;
                surf(-X, Y, Z);
                %plot3(-pnts(:,1), pnts(:,2), pnts(:,3));
            elseif (~xsym && ysym)
                hold on;
                surf(X, -Y, Z);
                %plot3(pnts(:,1), -pnts(:,2), pnts(:,3));
            elseif (xsym && ysym)
                hold on;
                surf(X, -Y, Z);
                surf(-X, Y, Z);
                surf(-X, -Y, Z);
%                 plot3(pnts(:,1), -pnts(:,2), pnts(:,3));
%                 plot3(-pnts(:,1), pnts(:,2), pnts(:,3));
%                 plot3(-pnts(:,1), -pnts(:,2), pnts(:,3));
            end
        end
        
        % Overloaded real operator
        function [panOut] = real(panIn)
            panOut = Panel(panIn.Vertices);
            
            if (~isempty(panIn.Value))
                panOut.Value = real(panIn.Value);
            end
        end
        
        % Overloaded imag operator
        function [panOut] = imag(panIn)
            panOut = Panel(panIn.Vertices);
            
            if (~isempty(panIn.Value))
                panOut.Value = imag(panIn.Value);
            end
        end
        
        % Overloaded abs operator
        function [panOut] = abs(panIn)
            panOut = Panel(panIn.Vertices);
            
            if (~isempty(panIn.Value))
                panOut.Value = abs(panIn.Value);
            end
        end
        
        % Overloaded angle operator
        function [panOut] = angle(panIn)
            panOut = Panel(panIn.Vertices);
            
            if (~isempty(panIn.Value))
                panOut.Value = angle(panIn.Value);
            end
        end
    end
    
    methods (Static, Access = private)
        function [pass] = checkVertices(verts)
            [n m] = size(verts);
            if (n ~= 4 || m ~= 3)
                error('The vertices must be a 4x3 matrix. i.e. 4 points by 3 coordinates (x,y,x)');
            end
            pass = true; 
        end
    end
    
    methods (Static, Access = public)
        function [norm, area, centroid] = Properties(verts)
            
            v1 = verts(2,:) - verts(1,:);
            v2 = verts(3,:) - verts(1,:);
            v3 = verts(4,:) - verts(1,:);
            
            % triangle 1
            n1 = cross(v1,v2);
            len = sqrt(n1(1)^2 + n1(2)^2 + n1(3)^2);
            ar1 = 0.5*len;
            n1 = n1./len;
            cen1 = mean(verts(1:3,:));
            
            % triangle 1
            n2 = cross(v1,v3);
            len = sqrt(n2(1)^2 + n2(2)^2 + n2(3)^2);
            ar2 = 0.5*len;
            n2 = n2./len;
            cen2 = mean(verts([1 3 4],:));
            
           
            
            area = ar1 + ar2;
            if (area == 0)
                norm = [0 0 0];
                centroid = mean(verts);
            else
                if (ar1 == 0)
                    norm = n2;
                elseif (ar2 == 0)
                    norm = n1;
                else
                    norm = mean([n1; n2]);
                end
                
                centroid = 1/area*(cen1*ar1 + cen2*ar2);
                
                % could check how coplanar the points are by comparing the
                % normals
            end
            
%             p = verts(3,:) - verts(1,:);
%             q = verts(4,:) - verts(2,:);
% 
%             n = cross(p,q); %[(p(2)*q(3)-p(3)*q(2)), -(p(1)*q(3)-p(3)*q(1)), (p(1)*q(2)-p(2)*q(1))];
%             len = sqrt(n(1)^2 + n(2)^2 + n(3)^2);
%             norm = n./len;           
            
           % area = 0.5*len;
        end
    end
    
    methods (Access = private)
        function [] = computeQuantities(pan)
            
            [norm, ar, cent] = Panel.Properties(pan.vertices);

            pan.normal = norm;
            pan.area = ar;
            
            pan.centroid = cent;
        end
    end
end