classdef Trapazoid
    
    properties (Access = private)
        panel;
    end
    
    properties (Dependent)
        Vertices;
        Area;
        Perimeter;
        Lengths;
        Centroid;
        Interia;
    end
    
    methods
        function [tr] = Trapazoid(verts)
            Trapazoid.checkVertices(verts);
            verts3 = zeros(4,3);
            verts3(:,1:2) = verts;
            tr.panel = Panel(verts3);
        end
        
        function [verts] = get.Vertices(tr)
            verts3 = tr.panel.Vertices;
            verts = verts3(:,1:2);
        end

        function [ar] = get.Area(tr)            
            ar = tr.panel.Area;
        end
        
        function [pr] = get.Perimeter(tr) 
            pr = sum(tr.Lengths);
        end
        
        function [lens] = get.Lengths(tr)  
            verts = tr.Vertices;  
            lens = zeros(4,1);
            for n = 1:3
                lens(n) = Trapazoid.distance(verts(n,:), verts(n+1,:));
            end
            lens(4) = Trapazoid.distance(verts(4,:), verts(1,:));
        end
        
        function [cen] = get.Centroid(tr)
            cen3 = tr.panel.Centroid;
            cen = [cen3(1) cen3(2)];
        end
        
        function [I] = get.Interia(tr)
            verts = tr.Vertices;            
            [~, ~, I] = Trapazoid.TrapAreaCen(verts);
        end
        
        function [xlim, ylim] = BoundBox(tr)
            verts = tr.Vertices;
            
            xlim(1) = min(verts(:,1));
            xlim(2) = max(verts(:,1));
            ylim(1) = min(verts(:,2));
            ylim(2) = max(verts(:,2));
        end
                
        function [] = Translate(tr, vector)
            if (length(vector) ~= 2)
                error('Translational vector must be 2x1 or 1x2');
            end
            
            vector3 = zeros(1,3);
            vector3(1:2) = vector(1:2);
            
            tr.panel.Translate(vector3);
        end
        
        function [trScale] = GetScaled(tr, scale)
            verts = tr.Vertices;
            cent = tr.Centroid;
            
            vertsCen = verts - ones(4,1)*cent;
            
            vertsScaleCen = scale*vertsCen;
            vertsScale = vertsScaleCen + ones(4,1)*scale*cent;
            
            trScale = Trapazoid(vertsScale);
        end
        
        function [xz] = Xz(tr, z)
            vs = tr.panel.Vertices;
            vs = vs(:,1:2);
            
            xz = zeros(1,2);
            
            for n = 1:2
                m = (vs(n,1)-vs(5-n,1))/(vs(n,2)-vs(5-n,2));
                xz(n) = m*(z - vs(5-n,2)) + vs(5-n,1);
            end
        end
        
        function [isIn] = IsInside(tr, pt)
            verts3 = tr.panel.Vertices;
            verts = verts3(:,1:2);
            verts(5,1:2) = verts(1,1:2);
            
            [in, on] = inpolygon(pt(1), pt(2), verts(:,1), verts(:,2));
            isIn = (in || on);
        end
    end
    
    methods (Static)
        function [tr] = MakeTrap(lenTop, lenBot, dLenF, height, draft)
            freeBrd = height - draft;
            x1 = [-lenTop/2, freeBrd];
            x2 = [lenTop/2, freeBrd];
            x4 = [-lenTop/2 + dLenF, -draft];
            x3 = [x4(1) + lenBot, -draft];
            
            tr = Trapazoid([x1; x2; x3; x4]);
        end
        
        function [A, C, Ig] = TrapAreaCen(verts)
            % Tri 1
            t1 = Triangle([verts(1,:); verts(2,:); verts(4,:)]);
            A1 = t1.Area;
            C1 = t1.Centroid;
            
            % Tri 1
            t2 = Triangle([verts(3,:); verts(4,:); verts(2,:)]);
            A2 = t2.Area;
            C2 = t2.Centroid;
            
            A = A1 + A2;
            C = 1/(A1 + A2)*(C1*A1 + C2*A2);
            

            Ig = 0;
            
            for n = 1:4
                if (n == 4)
                    tri = Triangle([verts(4,:); verts(1,:); C]);
                else
                    tri = Triangle([verts(n,:); verts(n+1,:); C]);
                end
                b = tri.Base;
                h = tri.Height;
                a = tri.Offset;
                % moment of interia formula for no oblique triangles
                I = 1/36*(b^3*h - b^2*h*a + b*h*a^2 + b*h^3);
                d = Trapazoid.distance(tri.Centroid, C);
                A = tri.Area;
                
                I = I + A*d^2;
                Ig = I + Ig;
            end
        end
    end
    
    methods (Static, Access = private)
        function [pass] = checkVertices(verts)
            [n, m] = size(verts);
            if (n ~= 4 || m ~= 2)
                error('The vertices must be a 4x2 matrix. i.e. 4 points by 2 coordinates (x,y)');
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