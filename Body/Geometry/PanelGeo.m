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
                        newPans(n).Value = pans.Panels(n).Value;
                        newPans(n).IsWet = pans.Panels(n).IsWet;
                        newPans(n).IsInterior = pans.Panels(n).IsInterior;
                        newPans(n).IsBody = pans.Panels(n).IsBody;
                    end
                    
                    geo.panels = newPans;
                    geo.xsym = pans.Xsymmetry;
                    geo.ysym = pans.Ysymmetry;
                else
                    geo.panels = pans;

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
        function [geo] = set.Values(geo, vals)
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
        
        function [] = Translate(geo, vector)
            N = length(geo.panels);
            for n = 1:N
                geo.panels(n).Translate(vector);
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
                
        function [] = plot(geo, varargin)
            geo.plotFuncs(@plot, varargin{:});
        end
        
        function [] = surf(geo, varargin)
            geo.plotFuncs(@surf, varargin{:});
        end
        
        % Overloaded plus operator
        function [geoOut] = plus(geoa, geob)
            
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
        function [] = plotFuncs(geo, func, varargin)
            N = geo.Count;
            
            opts = checkOptions({{'ShowSym'}, {'ShowNorm'}, {'OnlyWet'}}, varargin);
            showSym = opts(1);
            
            xsy = false;
            ysy = false;
            
            if (showSym)
                xsy = geo.xsym;
                ysy = geo.ysym;
            end
            
            if (~xsy && ~ysy)
                for n = 1:N
                    func(geo.panels(n), varargin{:});
                    hold on;
                end
            elseif (xsy && ~ysy)
                for n = 1:N
                    func(geo.panels(n), 'xsym', varargin{:});
                    hold on;
                end
            elseif (~xsy && ysy)
                for n = 1:N
                    func(geo.panels(n), 'ysym', varargin{:});
                    hold on;
                end
            elseif (xsy && ysy)
                for n = 1:N
                    func(geo.panels(n), 'xsym', 'ysym', varargin{:});
                    hold on;
                end
            end
        end
    end
end