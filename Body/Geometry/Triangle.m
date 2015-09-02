classdef Triangle
    
    properties (Access = private)
        verts;
        h;
    end
    
    properties (Dependent)
        Vertices;
        Area;
        Centroid;
        Base;
        Height;
        Offset;
    end
    
    methods
        function [tr] = Triangle(verts)
            Triangle.checkVertices(verts);
            
            tr.verts = verts;
            tr.h = Triangle.height(verts(1,:), verts(2,:), verts(3,:)); 
        end
        
        function [v] = get.Vertices(tr)
            v = tr.verts;
        end

        function [ar] = get.Area(tr)  
            v = tr.verts;
            v1 = v(2,:) - v(1,:);
            v2 = v(3,:) - v(1,:);
            
            v13 = [v1(1), v1(2), 0];
            v23 = [v2(1), v2(2), 0];
            
            n1 = cross(v13,v23);
            len = sqrt(n1(1)^2 + n1(2)^2 + n1(3)^2);
            ar = 0.5*len;
        end
        
        function [cen] = get.Centroid(tr)
            v = tr.verts;
            cen = 1/3*(v(1,:) + v(2,:) + v(3,:));
        end
        
        function [b] = get.Base(tr)
            v = tr.verts;
            b = Triangle.distance(v(1,:), v(2,:));
        end
        
        function [h] = get.Height(tr)
            v = tr.verts;
            h = Triangle.height(v(1,:), v(2,:), v(3,:)); 
        end
        
        function [off] = get.Offset(tr)
            v = tr.verts;
            l = Triangle.distance(v(1,:), v(3,:));
            off = sqrt(l^2 - tr.h^2);
        end
    end
    
    methods (Static, Access = private)
        function [pass] = checkVertices(verts)
            [n, m] = size(verts);
            if (n ~= 3 || m ~= 2)
                error('The vertices must be a 3x2 matrix. i.e. 3 points by 2 coordinates (x,y)');
            end
            pass = true; 
        end
        
        function [d] = distance(v1, v2)
            d = sqrt((v1(1)-v2(1))^2 + (v1(2)-v2(2))^2);
        end
        
        function [h] = height(l1, l2, p)
            dx = l2(1) - l1(1);
            dy = l2(2) - l1(2);            
            d = sqrt(dx^2 + dy^2);
            
            h = abs(dy*p(1) - dx*p(2) + l2(1)*l1(2) - l2(2)*l1(1))/d;
        end
    end
end