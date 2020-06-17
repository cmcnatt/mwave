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
        setNorm;
        normal;
        area;
        centroid;
        value;
        isWet;
        isIn;
        isBody;
    end
    
    properties (Dependent)
        Vertices;
        Normal;
        Area;
        Centroid;
        Value;
        IsWet;
        IsInterior;
        IsBody;
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
                pan.setNorm = false;
                pan.area = [];
                pan.centroid = [];
                pan.value = [];
                pan.isWet = true;
                pan.isIn = false;
                pan.isBody = true;
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
        function [] = set.Normal(pan, val)
            pan.setNorm = true;
            pan.normal = val;
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
        
        function [iw] = get.IsWet(pan)
            iw = pan.isWet;
        end
        function [pan] = set.IsWet(pan, iw)
            if (~isBool(iw))
                error('IsWet must be boolean');
            end
            
            pan.isWet = iw;
        end
        
        function [ii] = get.IsInterior(pan)
            ii = pan.isIn;
        end
        function [pan] = set.IsInterior(pan, ii)
            if (~isBool(ii))
                error('IsInterior must be boolean');
            end
            
            pan.isIn = ii;
        end
        
        function [ib] = get.IsBody(pan)
            ib = pan.isBody;
        end
        function [pan] = set.IsBody(pan, ib)
            if (~isBool(ib))
                error('IsBody must be boolean');
            end
            
            pan.isBody = ib;
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
        
        function [cent] = OffsetCentroid(pan, dist)
            cent = pan.Centroid;
            norm = pan.Normal;
            cent = cent + dist*norm;
        end
        
        function [] = Translate(pan, vector)
            if (length(vector) ~= 3)
                error('Translational vector must be 3x1 or 1x3');
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
            
            opts = checkOptions({{'xsym'}, {'ysym'}, {'ShowNorm'}, {'OnlyWet'}, {'ShowInt'}}, varargin);
            xsym = opts(1);
            ysym = opts(2);
            showN = opts(3);
            onlyWet = opts(4);
            showInt = opts(5);
            
            if (onlyWet)
                if (pan.isWet)
                    if (pan.isIn)
                        if (showInt)
                            plotThis = true;
                        else
                            plotThis =  false;
                        end
                    else
                        plotThis = true;
                    end
                else
                    plotThis = false;
                end
            else
                if (pan.isBody)
                    plotThis = true;
                else
                    plotThis = false;
                end
            end
                
            
%             plotThis = true;
%             if (onlyWet && ~pan.isWet)
%                 plotThis = false;                    
%             end
%             
%             if (~showInt && ~pan.IsBody && pan.IsInterior)
%                 plotThis = false;                    
%             end

            if (plotThis)
                pnts = zeros(5,3);
                for n = 1:3
                    pnts(:,n) = [pan.vertices(:,n); pan.vertices(1,n)];
                end
                plot3(pnts(:,1), pnts(:,2), pnts(:,3), 'Color', [0 0.4470 0.7410]);

                if (showN)
                    hold on;
                    cent = pan.Centroid;
                    norm = pan.Normal;
                    quiver3(cent(1), cent(2), cent(3), norm(1), norm(2), norm(3), 'Color', MColor.Black);
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
        end
        
        function [] = surf(pan, varargin)
            
            opts = checkOptions({{'xsym'}, {'ysym'}, {'ShowNorm'}, {'OnlyWet'}, {'ShowInt'}}, varargin);
            xsym = opts(1);
            ysym = opts(2);
            showN = opts(3);
            onlyWet = opts(4);
            showInt = opts(5);
            

            if (pan.IsInterior)
                if (showInt)
                    plotThis = true;
                else
                    plotThis = false;
                end
            else
                if (onlyWet)
                    if (pan.IsWet)
                        plotThis = true;
                    else
                        plotThis = false;
                    end
                else
                    plotThis = true;
                end
            end
            
%             if (pan.IsBody)
%                  if (onlyWet && ~pan.isWet && ~pan.IsInterior)
%                      plotThis = false;
%                  end
%             else
%                 if ~(showInt && pan.IsInterior)
%                     plotThis = false;
%                 end
%             end
            
            
            
            if (plotThis)
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
                error('The vertices must be a 4x3 matrix. i.e. 4 points by 3 coordinates (x,y,z)');
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
            if (len == 0)
                n1 = [0 0 0];
            else
                n1 = n1./len;
            end
            cen1 = sum(verts(1:3,:))./3;
            
            % triangle 2
            n2 = cross(v2,v3);
            len = sqrt(n2(1)^2 + n2(2)^2 + n2(3)^2);
            ar2 = 0.5*len;
            if (len == 0)
                n2 = [0 0 0];
            else
                n2 = n2./len;
            end
            cen2 = sum(verts([1 3 4],:))./3;

            area = ar1 + ar2;
            if (area == 0)
                norm = [0 0 0];
                centroid = sum(verts)./4;
            else
                if (ar1 == 0)
                    norm = n2;
                elseif (ar2 == 0)
                    norm = n1;
                else
                    norm = sum([n1; n2])./2;
                end
                
                centroid = 1/area*(cen1*ar1 + cen2*ar2);
                
                % could check how coplanar the points are by comparing the
                % normals
            end
        end
    end
    
    methods (Access = private)
        function [] = computeQuantities(pan)
            
            [norm, ar, cent] = Panel.Properties(pan.vertices);

            if ~pan.setNorm
                pan.normal = norm;
            end
            pan.area = ar;
            
            pan.centroid = cent;
        end
    end
end