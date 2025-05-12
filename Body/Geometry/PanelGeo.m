classdef PanelGeo < handle

    properties (Access = private)
        panels;
        xsym;
        ysym;
    end

    properties (Dependent)
        Panels;
        Count;
        Centroids;
        Normals;
        Areas;
        Values;
        IsWets;
        IsInteriors
        IsBodies;
        AvgPanelArea;
        Xsymmetry;
        Ysymmetry;
    end

    methods
        % Constructor
        function [geo] = PanelGeo(pans, varargin)
            if (nargin >= 1)
                if (isa(pans, 'PanelGeo'))
                    nPan = pans.Count;

                    for n = 1:nPan
                        newPans(n) = Panel(pans.Panels(n).Vertices);
                        newPans(n).Normal = pans.Panels(n).Normal;
                        newPans(n).Value = pans.Panels(n).Value;
                        newPans(n).IsWet = pans.Panels(n).IsWet;
                        newPans(n).IsInterior = pans.Panels(n).IsInterior;
                        newPans(n).IsBody = pans.Panels(n).IsBody;
                    end

                    geo.panels = newPans;
                    geo.xsym = pans.Xsymmetry;
                    geo.ysym = pans.Ysymmetry;
                else
                    N = length(pans);
                    for n = 1:N
                        newPans(n) = Panel(pans(n).Vertices);
                        newPans(n).Normal = pans(n).Normal;
                        newPans(n).Value = pans(n).Value;
                        newPans(n).IsWet = pans(n).IsWet;
                        newPans(n).IsInterior = pans(n).IsInterior;
                        newPans(n).IsBody = pans(n).IsBody;
                    end
                    geo.panels = newPans;

                    geo.xsym = false;
                    geo.ysym = false;
                    opts = checkOptions({'Xsym' 'Ysym'}, varargin);
                    geo.xsym = opts(1);
                    geo.ysym = opts(2);
                end
            end
        end

        function [pans] = get.Panels(geo)
            pans = geo.panels;
        end

        function [cnt] = get.Count(geo)
            cnt = length(geo.panels);
        end

        function [cents] = get.Centroids(geo)
            cnt = geo.Count;
            cents = zeros(cnt, 3);

            for n = 1:cnt
                cents(n,:) = geo.panels(n).Centroid;
            end
        end

        function [norms] = get.Normals(geo)
            cnt = geo.Count;
            norms = zeros(cnt, 3);

            for n = 1:cnt
                norms(n,:) = geo.panels(n).Normal;
            end
        end

        function [ars] = get.Areas(geo)
            cnt = geo.Count;
            ars = zeros(cnt, 1);

            for n = 1:cnt
                ars(n) = geo.panels(n).Area;
            end
        end

        function [iws] = get.IsWets(geo)
            cnt = geo.Count;
            iws = zeros(cnt, 1);

            for n = 1:cnt
                iws(n) = geo.panels(n).IsWet;
            end
        end

        function [iis] = get.IsInteriors(geo)
            cnt = geo.Count;
            iis = zeros(cnt, 1);

            for n = 1:cnt
                iis(n) = geo.panels(n).IsInterior;
            end
        end

        function [ibs] = get.IsBodies(geo)
            cnt = geo.Count;
            ibs = zeros(cnt, 1);

            for n = 1:cnt
                ibs(n) = geo.panels(n).IsBody;
            end
        end

        function [vals] = get.Values(geo)
            cnt = geo.Count;
            vals = zeros(cnt, 1);

            for n = 1:cnt
                valn = geo.panels(n).Value;
                if (isempty(valn))
                    vals(n) = NaN;
                else
                    vals(n) = valn;
                end
            end
        end
        function [] = set.Values(geo, vals)
            cnt = geo.Count;
            if (length(vals) ~= cnt)
                error('The number of values must be equal to the number of panels');
            end

            for n = 1:cnt
                geo.panels(n).Value = vals(n);
            end
        end

        function [area] = get.AvgPanelArea(geo)
            N = length(geo.panels);
            area = 0;
            for n = 1:N
                area = area + geo.panels(n).Area;
            end
            area = area/N;
        end

        function [xsy] = get.Xsymmetry(geo)
            xsy = geo.xsym;
        end

        function [ysy] = get.Ysymmetry(geo)
            ysy = geo.ysym;
        end

        function [] = set.Ysymmetry(geo,ysy)
            geo.ysym = ysy;
        end

        function [newGeo] = GetSubsetPanels(geo, incWet, incInt)
            pans = geo.panels;
            N = length(pans);
            inds = true(N, 1);
            if incWet
                inds = geo.IsWets == true;
            end
            if ~incInt
                inds = inds & (geo.IsInteriors == false);
            end

            newGeo = PanelGeo(pans(inds));
        end

        function [] = Translate(geo, vector)
            N = length(geo.panels);
            for n = 1:N
                geo.panels(n).Translate(vector);
            end
        end

        function [cents] = OffsetCentroids(geo, dist)
            cnt = geo.Count;
            cents = zeros(cnt, 3);

            for n = 1:cnt
                cents(n, :) = geo.panels(n).OffsetCentroid(dist);
            end
        end

        function [] = Rotate(geo, ax, angle, varargin)
            N = length(geo.panels);
            rotMat = createRotationMatrix(ax, angle);
            for n = 1:N
                geo.panels(n).Rotate('RotMat', rotMat, varargin{:});
            end
        end

        function [geo1, geo2] = Split(geo, planePnt, planeNorm)
            pt = planePnt;
            nrm = planeNorm;
            nrm = nrm./(nrm(1)^2 + nrm(2)^2 + nrm(3)^2);

            a = nrm(1);
            b = nrm(2);
            c = nrm(3);

            d = -dot(pt, nrm);

            N = length(geo.panels);

            nPos = 1;
            nNeg = 1;
            for m = 1:N
                panm = geo.panels(m);
                verts = panm.Vertices;
                allPos = true;
                for n = 1:4
                    v = verts(n,:);
                    D = a*v(1) + b*v(2) + c*v(3) + d;

                    if (D < 0)
                        allPos = false;
                    end
                end
                if (allPos)
                    pansPos(nPos) = Panel(verts);
                    pansPos(nPos).IsWet = panm.IsWet;
                    pansPos(nPos).IsBody = panm.IsBody;
                    pansPos(nPos).IsInterior = panm.IsInterior;
                    nPos = nPos + 1;
                else
                    pansNeg(nNeg) = Panel(verts);
                    pansNeg(nPos).IsWet = panm.IsWet;
                    pansNeg(nPos).IsBody = panm.IsBody;
                    pansNeg(nPos).IsInterior = panm.IsInterior;
                    nNeg = nNeg + 1;
                end
            end

            geo1 = PanelGeo(pansNeg);
            geo2 = PanelGeo(pansPos);
        end

        function [stl] = MakeStl(geo, cg, noInt)
            if nargin < 2
                cg = [0 0 0];
            end
            if nargin < 3
                noInt = true;
            end
            N0 = geo.Count;
            if noInt
                inclPan = ~geo.IsInteriors;
                N = N0 - sum(geo.IsInteriors);
            else
                inclPan = ones(N0, 1);
                N = N0;
            end
            verts0 = zeros(4*N, 3);
            faces0 = zeros(2*N, 3);

            im = 1;
            ifa = 1;
            for m = 1:N0
                if inclPan(m)
                    verts0(im:im+3,:) = geo.panels(m).Vertices;
                    faces0(ifa,:) = [im, im+1, im+2];
                    faces0(ifa+1,:) = [im+2, im+3, im];

                    im = im + 4;
                    ifa = ifa + 2;
                end
            end
            [verts, ~, ic] = unique(verts0, 'rows');

            faces = zeros(2*N, 3);
            for m = 1:2*N
                for n = 1:3
                    faces(m, n) = ic(faces0(m, n));
                end
            end

            stl = Stl6DOFGeo;

            stl.SetStl(verts, faces, cg);
        end

        function [] = WriteBpi(geo, fileLoc, name)
            fileName = [fileLoc '\' name '.bpi'];
            fileID = fopen(fileName, 'wt');

            ulen = 1;
            g = 9.806650;

            pans = geo.Panels;
            nPan = length(pans);
            fprintf(fileID, ['Model ' name ', created: ' date '\n']);
            fprintf(fileID, '%i\n', nPan);
            for n = 1:nPan
                fprintf(fileID, '\t%8.4f\t%8.4f\t%8.4f\n', verts(m,1), verts(m,2), verts(m,3));
            end

            fclose(fileID);
        end

        function [] = WriteGdf(geo, fileLoc, name, center, includeInt)
            if nargin < 4
                center = [];
                includeInt = true;
            end

            % make a copy
            geo0 = PanelGeo(geo);
            if ~isempty(center)
                geo0.Translate(-center); %translate to rotate about center of rotation.
            end

            filename = [fileLoc '\' name '.gdf'];
            fid = fopen(filename, 'wt');

            ulen = 1;
            g = 9.806650;

            Nwet = sum(geo0.IsWets);
            Nint = sum(geo0.IsInteriors);
            if (includeInt)
                count = Nwet;
            else
                count = Nwet - Nint;
            end

            fprintf(fid, ['Model ' name ', created: ' date '\n']);
            fprintf(fid, '%8.4f %8.4f\n', ulen, g);
            fprintf(fid, '%i %i \n', geo0.Xsymmetry, geo0.Ysymmetry);
            fprintf(fid, '%i\n', count);

            pans = geo0.Panels;

            for n = 1:geo0.Count
                pan = pans(n);

                panOk = true;
                if (~pan.IsWet)
                    panOk = false;
                end
                if (pan.IsInterior && ~includeInt)
                    panOk = false;
                end

                if (panOk)
                    verts = pan.Vertices;
                    for m = 1:4
                        fprintf(fid, '\t%8.4f\t%8.4f\t%8.4f\n', verts(m,1), verts(m,2), verts(m,3));
                    end
                end
            end

            fclose(fid);
        end

        function [faces, verts] = QuadMesh(geo, noInt, onlyWet)
            if nargin < 2
                noInt = true;
            end
            if nargin < 3
                onlyWet = false;
            end
            N0 = geo.Count;
            inclPan = geo.getIndices(noInt, onlyWet);
            N = sum(inclPan);
            verts0 = zeros(4*N, 3);
            faces0 = zeros(N, 4);

            im = 1;
            ifa = 1;
            for m = 1:N0
                if inclPan(m)
                    verts0(im:im+3,:) = geo.panels(m).Vertices;
                    faces0(ifa,:) = [im, im+1, im+2, im+3];

                    im = im + 4;
                    ifa = ifa + 1;
                end
            end
            [verts, ~, ic] = unique(verts0, 'rows');

            faces = zeros(N, 4);
            for m = 1:N
                for n = 1:4
                    faces(m, n) = ic(faces0(m, n));
                end
            end
        end

        function [] = plot(geo, varargin)
            geo.plotFuncs('plot', varargin{:});
        end

        function [] = surf(geo, varargin)
            geo.plotFuncs('surf', varargin{:});
        end

        function [geoOut] = plus(geoa, geob)
            % PLUS: overloade addition operator

            if (xor(geoa.Xsymmetry, geob.Xsymmetry) || xor(geoa.Ysymmetry, geob.Ysymmetry))
                error('PanelGeos must have the same symmetry to add');
            end

            nA = geoa.Count;
            nB = geob.Count;
            for n = 1:nA
                newPans(n) = Panel(geoa.Panels(n).Vertices);
            end

            for n = 1:nB
                newPans(n+nA) = Panel(geob.Panels(n).Vertices);
            end

            geoOut = PanelGeo(newPans);
        end

        function [geoOut] = times(arg1, arg2)
            % TIMES: overloaded times operator
            if isa(arg1, 'PanelGeo')
                geoIn = arg1;
                val = arg2;
            else
                geoIn = arg2;
                val = arg1;
            end

            geoOut = PanelGeo(geoIn);

            geoOut.Values = geoIn.Values.*val;
        end

        function [geoOut] = mtimes(arg1, arg2)
            % MTIMES: overloaded mtimes operator
            if isa(arg1, 'PanelGeo')
                geoIn = arg1;
                val = arg2;
            else
                geoIn = arg2;
                val = arg1;
            end

            geoOut = PanelGeo(geoIn);

            geoOut.Values = geoIn.Values*val;
        end

        function [geoOut] = rdivide(arg1, arg2)
            % RDIVIDE: overloaded rdivide operator
            if isa(arg1, 'PanelGeo')
                geoIn = arg1;
                val = arg2;
            else
                geoIn = arg2;
                val = arg1;
            end

            geoOut = geoIn.*(1./val);
        end

        % Overloaded real operator
        function [geoOut] = real(geoIn)
            geoOut = PanelGeo(geoIn);

            geoOut.Values = real(geoIn.Values);
        end

        % Overloaded imag operator
        function [geoOut] = imag(geoIn)
            geoOut = PanelGeo(geoIn);

            geoOut.Values = imag(geoIn.Values);
        end

        % Overloaded abs operator
        function [geoOut] = abs(geoIn)
            geoOut = PanelGeo(geoIn);

            geoOut.Values = abs(geoIn.Values);
        end

        % Overloaded angle operator
        function [geoOut] = angle(geoIn)
            geoOut = PanelGeo(geoIn);

            geoOut.Values = anlge(geoIn.Values);
        end
    end

    methods (Access = private)
        function inclPan = getIndices(obj, noInt, onlyWet)
            inclPan = true(obj.Count, 1);
            if noInt
                inclPan = inclPan & ~obj.IsInteriors;
            end
            if onlyWet
                inclPan = inclPan & obj.IsWets;
            end
        end

        function plotFuncs(obj, func, varargin)

            [opts, args] = checkOptions({{'ShowNorm'}, {'OnlyWet'}, {'NoInt'}, {'color', 1}, {'alpha', 1}, {'edge', 1}, {'line', 1}, {'water'}}, varargin);
            showNorm = opts(1);
            onlyWet = opts(2);
            noInt = opts(3);
            color = [];
            if opts(4)
                color = args{4};
            end
            alpha = 1;
            if opts(5)
                alpha = args{5};
            end
            edge = [];
            if opts(6)
                edge = args{6};
            end
            line = [];
            if opts(7)
                line = args{7};
            end
            showWater = opts(8);


            xsy = false;
            ysy = false;

            [mesh, verts] = obj.QuadMesh(noInt, onlyWet);

            vals = obj.Values;

            args = {'faces', mesh, 'vertices', verts,...
                'facelighting', 'none', 'facecolor', 'flat', 'edgelighting', 'flat',...
                'parent',gca};

            if strcmp(func, 'plot')
                colV = MColor.Blue;
                lineW = 0.5;
                if ~isempty(edge)
                    colV = edge;
                end
                if ~isempty(line)
                    lineW = line;
                end
                args{length(args) + 1} = 'edgecolor';
                args{length(args) + 1} = colV;
                args{length(args) + 1} = 'facecolor';
                args{length(args) + 1} = 'none';
                args{length(args) +1 } = 'linewidth';
                args{length(args) +1 } = lineW;
            elseif strcmp(func, 'surf')
                colV = MColor.Black;
                args{length(args) + 1} = 'edgecolor';
                args{length(args) + 1} = colV;
                args{length(args) + 1} = 'facealpha';
                args{length(args) + 1} = alpha;

                if ~isempty(color)
                    args{length(args) + 1} = 'facecolor';
                    args{length(args) + 1} = color;
                else
                    if isnan(vals(1))
                        cents = obj.Centroids;
                        colF = cents(:,3);
                    else
                        colF = vals;
                    end
                    args{length(args) + 1} = 'facevertexcdata';
                    args{length(args) + 1} = colF;
                    args{length(args) + 1} = 'cdata';
                    args{length(args) + 1} = colF;
                end

                if ~isempty(edge)
                    args{length(args) + 1} = 'edgecolor';
                    args{length(args) + 1} = edge;
                end
            end

            patch(args{:});

            if (showNorm)
                hold on;
                ipan = obj.getIndices(noInt, onlyWet);
                cent = obj.Centroids;
                norm = obj.Normals;
                area = obj.Areas;
                quiver3(cent(ipan,1), cent(ipan,2), cent(ipan,3), norm(ipan,1), norm(ipan,2), norm(ipan,3),'color', MColor.Black);
                %                 h1 = quiver3(cent(ipan,1), cent(ipan,2), cent(ipan,3), norm(ipan,1), norm(ipan,2), norm(ipan,3),'color', MColor.Black,'linewidth',1.5);
                %                 set(h1,'AutoScale','on', 'AutoScaleFactor', 1.5)
            end

            if showWater
                % Find maximum extents of geometry
                xExtent(1) = min(obj.Centroids(:,1));
                xExtent(2) = max(obj.Centroids(:,1));
                yExtent(1) = min(obj.Centroids(:,2));
                yExtent(2) = max(obj.Centroids(:,2));

                del = 0.4*max([-xExtent(1) xExtent(2) -yExtent(1) yExtent(2)]);

                water = PanelGeo(...
                    Panel(...
                    [xExtent(1)-del yExtent(1)-del 0;
                    xExtent(2)+del yExtent(1)-del 0;
                    xExtent(2)+del yExtent(2)+del 0;
                    xExtent(1)-del yExtent(2)+del 0]));
                [meshWater, vertsWater] = water.QuadMesh(noInt, onlyWet);
                argsWater = args; % Reuse args from WEC mesh, but change some settings
                argsWater{2} = meshWater; % faces
                argsWater{4} = vertsWater; % vertices
                argsWater{14} = 'none'; % turn off edges
                argsWater{16} = 0.5; % Add transparency
                argsWater{18} = MColor.Blue; % Set colour
                patch(argsWater{:});
                %                 surf(water, 'color', MColor.Blue, 'alpha', 0.5);
            end

        end
    end
end