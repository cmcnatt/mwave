classdef FloatingCylinder < FloatingBody
    
    properties (Access = protected)
        radius;
        draft;
        height;
    end

    properties (Dependent)
        Radius;
        Height;
        Draft;
    end
    
    methods
        % Constructor 
        function [fb] = FloatingCylinder(rho, radius, height, draft, varargin)
            % Floating cylinder of uniform density that rotates about its
            % cg
            %
            % rho - fluid density
            % radius - radius of the cylinder
            % height - height of the cylinder from draft to extension above
            % water
            % draft - draft of subergence
            % optional arguments - number or panels in theta (Ntheta), r
            % (Nr), and z (Nz)
            
            fb = fb@FloatingBody();
            zcg0 = -draft+height/2;
            fb.cg = [0 0 zcg0];

            fb.m = FloatingCylinder.MassMatrix(rho, radius, height, draft);
            fb.position = [0 0 0];
            
            fb.radius = radius;
            fb.height = height;
            fb.draft = draft;
            
            if (~isempty(varargin))
                if (length(varargin) < 3)
                    error('There must be three optional inputs to FloatingCylinder: Ntheta, Nr, Nz');
                end
                Ntheta = varargin{1};
                Nr = varargin{2};
                Nz = varargin{3};
                
                opts = checkOptions({{'UseSym'}, {'NoInt'}}, varargin);
                optsIn = {};
                n = 0;
                if (opts(1))
                    n = n + 1;
                    optsIn{n} = 'Quarter';
                end
                if (opts(2))
                    n = n + 1;
                    optsIn{n} = 'NoInt';
                    fb.iSurfPan = 0;
                else
                    fb.iSurfPan = 1;
                end
                
                panGeo = makePanel_cylinder(fb.radius, fb.draft, fb.height, Ntheta, Nr, Nz, optsIn{:});
                fb.panelGeo = panGeo;
                fb.iLowHi = 0;
            end
            
            pts = CirContSurf.Compute(fb.cg, radius, 30);
            pts = [pts; radius 0 0];
            fb.wpSec = pts(:, 1:2);
        end
    end
    
    methods
        function [r] = get.Radius(fb)
            r = fb.radius;
        end
        
        function [h] = get.Height(fb)
            h = fb.height;
        end
        
        function [d] = get.Draft(fb)
            d = fb.draft;
        end
        
        function [] = MakePanelGeometry(fb, Ntheta, Nr, Nz, varargin)
            opts = checkOptions({'UseSym'}, varargin);
            
            if (opts(1))
                panGeo = makePanel_cylinder(fb.radius, fb.draft, fb.draft, Ntheta, Nr, Nz, 'Quarter');
            else
                panGeo = makePanel_cylinder(fb.radius, fb.draft, fb.draft, Ntheta, Nr, Nz);
            end
            %panGeo.Translate(-fb.position);
            fb.panelGeo = panGeo;
            fb.iLowHi = 0;
        end
    end
    
    methods (Access = protected)
        function [] = onModifyCg(fb, cg)
            mass = fb.m(1);
            wetVol = pi*fb.radius^2*fb.draft;
            rho = mass/wetVol;
            zcg = cg(3);

            fb.m = FloatingCylinder.MassMatrix(rho, fb.radius, fb.height, fb.draft, zcg);
        end
    end
    
    methods (Static)
        function [M] = MassMatrix(rho, radius, height, draft, varargin)
            if (~isempty(varargin))
                % zcg is the vertical position in body coordinates, 
                % where z=0 is the calm water free surface
                zcg0 = -draft+height/2;
                zcg = varargin{1};
                delzcg = zcg - zcg0;
            else
                delzcg = 0;
            end
            % rho is the fluid density
            % body mass equals buoyant force
            if (draft > height)
                % assume neutrally buoyant
                draft = height;
            end
            wetVol = pi*radius^2*draft;
            m = rho*wetVol;

            Izz = m*radius^2/2;
            Ixx = 1/12*m*(3*radius^2 + height^2);
            Iyy = Ixx;
            
            if (delzcg ~= 0)
                Ixx = Ixx + m*delzcg^2;
                Iyy = Iyy + m*delzcg^2;
            end

            M = zeros(6,6);

            M(1,1) = m;
            M(2,2) = m;
            M(3,3) = m;
            M(4,4) = Ixx;
            M(5,5) = Iyy;
            M(6,6) = Izz;
        end
    end
end