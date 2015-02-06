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
classdef WaveField < IWaveField & matlab.mixin.Heterogeneous
    % A wave field provides information (pressure, elevation, velocity,
    % significant wave height, spectrums) at discrete points in space
    % specified on either a grid (array) or at arbitrary (x,y,z) points.  
    %
    % For now, the array option only provides information at the surface (z
    % = 0) for an array of X and Y
    % 
    % The generic wave field class contains different wave fields for
    % different wave periods, but not different directions.  Different wave
    % direction of the same period are coherent and so in some cases should
    % be combined.  To handle wave field with seperate directional
    % components, use the DirWaveField class, which contains seperate
    % WaveField objects for each direction.
    %
    % The inputs to the base WaveField class are 
    %   - rho: the water density
    %   - g: the gravitational constant
    %   - h: the water depth
    %   - t: the vector of periods
    %   - p: a matrix of the pressure values (For an array, it must be of
    %   size [Nt x Nx x Ny] where Nt is the number of periods, Nx is the
    %   number of points in the x-direction and Ny is the number of points
    %   in the y-direction. For a list of points, it must be [Nt x Np]
    %   where Np is the number of points.
    %   - v: a matris of the velocity values (For an array it must be of
    %   size [Nt x 3 x Nx x Ny], and for a list of points, it must be [Nt x
    %   3 x Np]
    %   - isarray: indicates whether the wave field is an array or a list
    %   of points
    %   - If isarray, then also must contain a mesh grid of X and Y
    %   - If not is array, then must contain a [Np x 3] vector of {x,y,z}
    %   points
    
    properties (Access = private)
        p;              % Nt x Ny x Nx (Nt x Npoints)
        vel;            % Nt x 3 x Ny x Nx (Nt x Nbeta x 3 x Npoints)
    end

    methods
        
        % Constructor
        function [wf] = WaveField(rho, g, h, t, p, vel, isarray, varargin)
            if (nargin > 0)
                wf.isarray = isarray;
                wf.rho = rho;
                wf.g = g;
                [n m] = size(t);
                if (n ~= 1 && m ~= 1)
                    error('The periods must be an Nx1 or a 1xN vector');
                end
                if (n == 1)
                    t = t';
                end
                wf.t = t;
                wf.nT = length(t);
                wf.h = h;

                if (isarray)
                    % TODO: arrays only handle grids of X and Y, may want
                    % to include Z
                    X = varargin{1};
                    Y = varargin{2};
                    [nye, nxe] = size(X);
                    wf.x = X;
                    wf.y = Y;
                    
                    if (iscell(p))
                        nt = size(p,1);
                        [ny, nx] = size(p{1});
                        if (nt ~= wf.nT || ny ~= nye || nx ~= nxe)
                            error('Pressure input is not of the expected size');
                        end
                                       
                        wf.p = p;
                        
                        if (~isempty(vel))
                            if (~iscell(vel))
                                error('The pressure and velocity inputs to the wave field must be in the same format, cell or not cell');
                            end

                            [nt, nc] = size(vel);
                            [ny, nx] = size(vel{1,1});
                            if (nt ~= wf.nT || ny ~= nye || nx ~= nxe || nc ~= 3)
                                error('Velocity input is not of the expected size');
                            end

                            wf.hasvel = true;
                            wf.vel = vel;
                        else
                            wf.hasvel = false;
                        end
                    else
                        [nt, ny, nx] = size(p);
                        if (nt ~= wf.nT || ny ~= nye || nx ~= nxe)
                            error('Pressure input is not of the expected size');
                        end
                        
                        if (~isempty(vel))
                            wf.hasvel = true;
                            if (iscell(vel))
                                error('The pressure and velocity inputs to the wave field must be in the same format, cell or not cell');
                            end
                            
                            [nt, nc, ny, nx] = size(vel);
                            if (nt ~= wf.nT || ny ~= nye || nx ~= nxe || nc ~= 3)
                                error('Velocity input is not of the expected size');
                            end
                            
                            wf.vel = cell(wf.nT,3);
                        else
                            wf.hasvel = false;
                        end

                        wf.p = cell(wf.nT,1);
                        
                        for n = 1:wf.nT
                            wf.p{n} = squeeze(p(n,:,:));
                            if (wf.hasvel)
                                for m = 1:3
                                    wf.vel{n,m} = squeeze(vel(n,m,:,:));
                                end
                            end
                        end
                    end
                else
                    pts = varargin{1};
                    npe = size(pts, 1);
                    wf.points = pts;
                    if (iscell(p))
                        nt = size(p,1);
                        np = length(p{1});
                        if (nt ~= wf.nT || np ~= npe)
                            error('Pressure input is not of the expected size');
                        end
                        
                        wf.p = p;
                        
                        if (~isempty(vel))
                            if (~iscell(vel))
                                error('The pressure and velocity inputs to the wave field must be in the same format, cell or not cell');
                            end

                            [nt, nc] = size(vel);
                            np = length(vel{1,1});
                            if (nt ~= wf.nT || np ~= npe || nc ~= 3)
                                error('Velocity input is not of the expected size');
                            end
                                                    
                            wf.vel = vel;
                            wf.hasvel = true;
                        else
                            wf.hasvel = false;
                        end
                    else
                        [nt, np] = size(p);
                        if (nt ~= wf.nT || np ~= npe)
                            error('Pressure input is not of the expected size');
                        end
                        
                        wf.p = cell(wf.nT,1);
                        
                        if (~isempty(vel))
                            if (iscell(vel))
                                error('The pressure and velocity inputs to the wave field must be in the same format, cell or not cell');
                            end

                            if (ndims(vel) == 2)
                                [nc, nc] = size(vel);
                            else
                                [nt, nc, np] = size(vel);
                            end
                            if (nt ~= wf.nT || np ~= npe || nc ~= 3)
                                error('Velocity input is not of the expected size');
                            end
                            
                            wf.vel = cell(wf.nT,3);
                            wf.hasvel = true;
                        else
                            wf.hasvel = false;
                        end

                        for n = 1:wf.nT
                            wf.p{n} = squeeze(p(n,:));
                            if (wf.hasvel)
                                for m = 1:3
                                    wf.vel{n,m} = squeeze(vel(n,m,:)).';
                                end
                            end
                        end
                    end
                end
            end
        end
                
        % Pressure, Velocity, Elevation, Spectrums, SigWaveHeight
                                
        % Pressure
        function [pr] = Pressure(wf, varargin)
            % complex pressure values at wave field points.  Optional input
            % of 'Surface' available which returns only pressure values at
            % z = 0
            surf = 0;
            for n = 1:length(varargin)
                if (strcmp(varargin{n}, 'Surface'))
                    surf = 1;
                end
            end
            
            if (~wf.isarray && surf)
                cnt = 0;
                for n = 1:length(wf.points)
                    if (wf.points(n,3) == 0)
                        cnt = cnt + 1;
                        indSurf(cnt) = n;
                    end
                end
                
                pr = cell(wf.nT,1);
                for n = 1:wf.nT
                    pr{n} = wf.p{n}(indSurf);
                end
            else
                pr = wf.p;
            end
        end
                
        % Elevation
        function [eta] = Elevation(wf, varargin)
            % complex wave elevation values at wave field points on surface
            % (z = 0)
            pr = wf.Pressure('Surface');
            
            eta = cellfun(@(x) x./(wf.rho*wf.g), pr, 'UniformOutput', false);
        end
        
        % Velocity
        function [v] = Velocity(wf, varargin)
            % complex velocity values at wave field points
            if (wf.hasvel)
                v = wf.vel;
            else
                error('Wave field does not contain velocity');
            end
        end
        
        % SigWaveHeight
        function [hs] = SigWaveHeight(wf, varargin)
            % spectral significant wave height at surface field points
            eta = wf.Elevation;
            
            if (wf.isarray)
                hs = zeros(size(wf.x));
            else
                hs = zeros(size(wf.points,1));
            end
            
            for n = 1:wf.nT
                hs(:,:) = hs(:,:) + 0.5*abs(eta{n}).^2;
            end
                
            hs = 4*sqrt(hs);
        end
        
        function [specs, actPoints] = Spectra(wf, varargin)
            % spectra at either all points in the wave field or specified
            % points (optional argument 'Points', points (Nx2)
            
            usePts = 0;
            pts = [];
            for n = 1:length(varargin)
                if (strcmp(varargin{n}, 'Points'))
                    usePts = 1;
                    pts = varargin{n+1};
                end
            end
            
            eta = wf.Elevation;

            f = 1./wf.t;
            df = computeDelta(f);
            
            if (usePts)
                if (~wf.isarray)
                    error('Point option for Spectrum is only valid in an array wave field.');
                end
                
                [npts col] = size(pts);
                if (col ~= 2)
                    error('points must be an Nx2 array');
                end

                actPoints = NaN(npts, 2);
                
                specs(npts,1) = WaveSpectrum;

                for n = 1:npts
                    xi = pts(n,1);
                    yi = pts(n,2);
                    [xo ix yo iy] = wf.FindClosestPoint(xi, yi);
                    actPoints(n,:) = [xo yo];

                    a_ = zeros(wf.nT,1);
                    for m = 1:wf.nT
                        a_(m) = abs(eta{m}(iy, ix));
                    end

                    S = (0.5*a_.^2./df)';

                    specs(n) = WaveSpectrum(S, f);
                end
            else
                if (wf.isarray)
                    npts = wf.x(1,:);
                    yv = wf.y(:,1);

                    specs(length(yv), length(npts)) = WaveSpectrum;

                    for n = 1:length(yv)
                        for m = 1:length(npts)

                            a_ = zeros(wf.nT,1);
                            for o = 1:wf.nT
                                a_(o) = abs(eta{o});
                            end

                            S = (0.5*a_.^2./df)';

                            specs(n, m) = WaveSpectrum(S, f);
                        end
                    end
                else
                    npts = length(eta{1});
                    specs(npts,1) = WaveSpectrum;

                    for n = 1:npts
                        a_ = zeros(wf.nT,1);
                        for m = 1:wf.nT
                            a_(m) = abs(eta{m}(n));
                        end

                        S = (0.5*a_.^2./df)';

                        specs(n) = WaveSpectrum(S, f);
                    end
                end
                
                actPoints = [];
            end
        end
        
        % Flux
        
        function [flux] = EnergyFlux(wf, surf, varargin)
            if (~wf.hasvel)
                error('Cannot compute flux. Wave field does not contain velocity components.');
            end

            if (~isa(surf, 'ControlSurface'))
                error('surf must be a ControlSurface.');
            end
            
            total = true;
            
            for n = 1:length(varargin)
                switch varargin{n}
                    case 'Total'
                        total = true;
                    case 'PerPoint'
                        total = false;
                end
            end
            
            pts = surf.Points';
            nrms = surf.Norms';
            ars = surf.Areas;

            if (size(ars, 2) == 1)
                ars = ars';
            end

            npts = size(pts,2);
            
            wcs = PlaneWaves(ones(size(wf.t)), wf.t, zeros(size(wf.t)), wf.h);
                            
            Flux = zeros(wf.nT, npts);
            
            if (wf.isarray)
                if (~surf.Is2dXSection)
                    error('Because the WaveField array points are on the z = 0 surface, the ControlSurface must be a 2D cross section');
                end
                
                X = wf.x;
                Y = wf.y;
                
                Pr = wf.Pressure;
                Velo = wf.Velocity;

                P = cell(wf.nT, 1);
                Vel = cell(wf.nT, 3);
                
                for n = 1:wf.nT
                    P{n} = interp2(X, Y, Pr{n}, pts(1,:), pts(2,:));
                    Vel{n,1} = interp2(X, Y, Velo{n,1}, pts(1,:), pts(2,:));
                    Vel{n,2} = interp2(X, Y, Velo{n,2}, pts(1,:), pts(2,:));
                    Vel{n,3} = interp2(X, Y, Velo{n,3}, pts(1,:), pts(2,:));
                end
            else
%                 if any(surf.Points ~= wf.points)
%                     error('For non-array wave fields, the points in the wave field must match the points on the surface');
%                 end
                
                P = wf.Pressure;
                Vel = wf.Velocity;
            end
            
            if (surf.Is2dXSection)
                coeffs = 1/4*1/IWaves.G*[wcs.C].*[wcs.Cg];         
                %Coeffs = coeffs*ars;
            else
                coeffs = 1/4*ones(wf.nT,1);
                %Coeffs = coeffs*ars;
            end
            
            for n = 1:wf.nT
                vdotn = Vel{n,1}.*nrms(1,:) + Vel{n,2}.*nrms(2,:) + Vel{n,3}.*nrms(3,:);

                Flux(n,:) = coeffs(n)*(conj(P{n}).*vdotn + P{n}.*conj(vdotn));
            end
            
            if (total)
                %Flux = sum(Flux, 2);
                Flux = Flux*ars.';
            end
            
            flux = cell(wf.nT, 1);
            for n = 1:wf.nT
                flux{n} = Flux(n,:);
            end
        end
        
        function [] = RemoveGeometries(wf, geos, wfInOrOut)
            
            if (~iscell(geos))
                geos = {geos};
            end
            
            nGeo = length(geos);
            
            if (strcmp(wfInOrOut, 'In'))
                wfIn = true;
            elseif (strcmp(wfInOrOut, 'Out'))
                wfIn = false;
            else
                error('The wfInOrOut value must be either ''In'' or ''Out''.');
            end
            
            for m = 1:nGeo
                geon = geos{m};
                
                [M, N] = size(geon);
                
                if (M < 4)
                    error('The geometry must have at least 4 (x,y) points');
                end
                
                if (N ~= 2)
                    error('The geometry is defined by an Nx2 matrix of N (x,y) points');
                end
                
                if all(geon(1,:) ~= geon(end,:))
                    error('The first and last point of the geometry must be same to create a closed contour');
                end
                
                if (wf.isarray)
                    xq = wf.x;
                    yq = wf.y;
                else
                    xq = wf.points(:,1);
                    yq = wf.points(:,2);
                end
                
                [in, on] = inpolygon(xq, yq, geon(:,1),geon(:,2));
                
                if (wfIn)
                    in = ~in;
                end
                
                for n = 1:wf.nT
                    pm = wf.p{m};
                    pm(in) = NaN;
                    wf.p{m} = pm;
                    
                    if (wf.hasvel)
                        velm = wf.vel{m};
                        velm(in,:) = NaN;
                        wf.vel{m} = velm;
                    end
                end
            end
        end
                        
        % Overloaded operators
        
        % Overloaded equality operator
        function [areEq] = eq(wfa, wfb)
            areEq = 1;
            
            if (~isa(wfa, 'IWaveField') || ~isa(wfb, 'IWaveField'))
                error('Each argument must be a IWaveField');
            end
            
            if (isa(wfb, 'WaveFieldCollection'))
                if (wfb == wfa)
                    areEq = 1;
                    return;
                else
                    areEq = 0;
                    return;
                end
            end
            
            % rho
            if (wfa.Rho ~= wfb.Rho)
                areEq = 0;
                return;
            end
            
            % array
            if (wfa.IsArray ~= wfb.IsArray)
                areEq = 0;
                return;
            end
            
            % points
            if (wfa.IsArray)
                [Xa, Ya] = wfa.FieldPoints;
                [Xb, Yb] = wfa.FieldPoints;
                
                if (any(any(Xa ~= Xb)) || any(any(Ya ~= Yb)))
                    areEq = 0;
                    return;
                end
            else
                pa = wfa.FieldPoints;
                pb = wfa.FieldPoints;
                
                if (any(any(pa ~= pb)))
                    areEq = 0;
                    return
                end
            end
            
            % frequncies
            if any(wfa.T ~= wfb.T)
                areEq = 0;
                return;
            end
            
            % depth
            if (wfa.H ~= wfb.H)
                areEq = 0;
                return;
            end
        end
        
        % Overloaded inequality operator
        function [areNE] = ne(wfa, wfb)
            areNE = 1;
            if (wfa == wfb)
                areNE = 0;
            end
        end
        
        % Overloaded plus operator
        function [wfout] = plus(wfa, wfb)
            if (wfa ~= wfb)
               error('Wave fields do not have the same density, points, or periods.  Cannot add.');
            end
            
            if (isa(wfb, 'WaveFieldCollection'))
                wfout = wfb + wfa;
                return;
            end
            
            plus = @(x,y) x + y;

            pa = wfa.Pressure;
            pb = wfb.Pressure;
            pout = cellfun(plus, pa, pb, 'UniformOutput', false);

            if (wfa.hasvel && wfb.hasvel)
                vela = wfa.Velocity;
                velb = wfb.Velocity;
                velout = cellfun(plus, vela, velb, 'UniformOutput', false);
            else
                velout = [];
            end

            if (wfa.IsArray)
                [X Y] = wfa.FieldPoints;
                wfout = WaveField(wfa.Rho, wfa.G, wfa.H, wfa.T, pout, velout, 1, X, Y);
            else
                pts = wfa.FieldPoints;
                wfout = WaveField(wfa.Rho, wfa.G, wfa.H, wfa.T, pout, velout, 0, pts);
           end
        end
                
        % Overloaded unary minus operator
        function [wfout] = uminus(wfin)
            neg = @(x) -x;
            pout = cellfun(neg, wfin.Pressure, 'UniformOutput', false);
            if (wfin.hasvel)
                velout = cellfun(neg, wfin.Velocity, 'UniformOutput', false);
            else
                velout = [];
            end
            
            if (wfin.IsArray)
               [X Y] = wfin.FieldPoints;
               wfout = WaveField(wfin.Rho, wfin.G, wfin.H, wfin.T, pout, velout, 1, X, Y);
            else
                pts = wfin.FieldPoints;
               wfout = WaveField(wfin.Rho, wfin.G, wfin.H, wfin.T, pout, velout, 0, pts);
           end
        end
        
        % Overloaded unary minus operator
        function [wfout] = minus(wfa, wfb)
            wfb = -wfb;
            wfout = wfa + wfb;
        end
        
        % Overloaded matrix multiplication operator
        function [wfout] = mtimes(a, b)
            if (isnumeric(a) && isa(b, 'WaveField'))
                num = a;
                wfin = b;
            elseif (isnumeric(b) && isa(a, 'WaveField'))
                num = b;
                wfin = a;
            else
                error('Multiplication only defined between a number and a WaveField.')
            end
            
            pin = wfin.Pressure;
            if (wfin.hasvel)
                velin = wfin.Velocity;
            else
                velout = [];
            end
            
            if (length(num) > 1)
                if (length(num) ~= wfin.nT)
                    error('The size of the vector must be the same as the number of periods in the wave field');
                end
                
                pout = cell(wfin.nT,1);
                if (wfin.hasvel)
                    velout = cell(wfin.nT,3);
                end
                
                for n = 1:wfin.nT
                    pout{n} = num(n)*pin{n};
                    pin{n} = [];
                    if (wfin.hasvel)
                        for m = 1:3
                            velout{n,m} = num(n)*velin{n,m};
                        end
                    end
                end
            else
                mult = @(x) num*x;
                    
                pout = cellfun(mult, pin, 'UniformOutput', false);
                if (wfin.hasvel)
                    velout = cellfun(mult, velin, 'UniformOutput', false);
                end
            end
            
            if (wfin.IsArray)
               [X, Y] = wfin.FieldPoints;
               wfout = WaveField(wfin.Rho, wfin.G, wfin.H, wfin.T, pout, velout, 1, X, Y);
            else
               pts = wfin.FieldPoints;
               wfout = WaveField(wfin.Rho, wfin.G, wfin.H, wfin.T, pout, velout, 0, pts);
           end
        end
        
        % Overloaded multiplication operator
        function [wfout] = times(a, b)
            wfout = a*b;
        end
    end    
end