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
classdef BodySurfWaveField < FBWaveField
    
    properties (Access = private)
        bodyGeo;
        bodyInds;
        motionFuncs;
        motBodyInds;
    end
    
    properties (Dependent)
        BodyGeo;
        BodyInds;
        BodyCount;
        MotionFuncs;
        MotionBodyInds;
    end
    
    methods
        function [wf] = BodySurfWaveField(body, iwave, swave, varargin)
            
            if (~isa(body, 'PanelGeo'))
                error('Body must be a PanelGeo');
            end
            if (iwave.IsArray)
                error('BodySurfWaveField must be a collection of points');
            end
            
            bodyPnts = body.Centroids;
            wfPnts = iwave.FieldPoints;
            
            Npnts = size(bodyPnts, 1);
            
            if (Npnts ~= size(wfPnts, 1))
                error('The body points and wave field points must be the same');
            end
            
            for n = 1:Npnts
                if any(abs(bodyPnts(n,:) - wfPnts(n,:)) > 1e-6*[1 1 1])
   %                  error('The body points and wave field points must be the same');
                end
            end
            
            wf = wf@FBWaveField(iwave, swave, varargin{:});
            
            wf.bodyGeo = body;
        end
        
        function [bod] = get.BodyGeo(wf)
            bod = wf.bodyGeo;
        end
        
        function [val] = get.BodyInds(wf)
            val = wf.bodyInds;
        end
        function [] = set.BodyInds(wf, val)
            if length(val) ~= wf.bodyGeo.Count
                error('BodyInds vector must be same length as number of panels');
            end
            wf.bodyInds = val;
        end
        
        function [val] = get.BodyCount(wf)
            val = length(unique(wf.bodyInds));
        end
        
        function [val] = get.MotionFuncs(wf)
            val = wf.motionFuncs;
        end
        function [] = set.MotionFuncs(wf, val)
            if ~isa(val, 'IMotionFunc')
                error('MotionFuncs must be an IMotionFunc');
            end
            wf.motionFuncs = val;
        end
        
        function [val] = get.MotionBodyInds(wf)
            val = wf.motBodyInds;
        end
        function [] = set.MotionBodyInds(wf, val)
            if ~isempty(wf.motionFuncs)
                if length(val) ~= length(wf.motionFuncs)
                    error('MotionBodyInds vector must be same length as number of MotionFuncs');
                end
            end
            wf.motBodyInds = val;
        end
                
        % Pressure
        function [pgeo] = Pressure(wf, type)
            
            rho = wf.Rho;
            g = wf.G;
            points = wf.FieldPoints;
            [Np, ~] = size(points);
            
            if strcmpi(type, 'Hydrostatic')
                z = points(:,3);
                z(z > 0) = 0;
                p = {-rho*g*z};
            elseif strcmpi(type, 'MotionHydrostatic')
                if isempty(wf.motionFuncs)
                    error('WaveField.MotionFuncs empty. Cannot compute MotionHydrostatic');
                end
                
                motions = wf.BodyMotions;
                funcs = wf.motionFuncs;
                
                bInds = wf.bodyInds;
                mInds = wf.motBodyInds;
                
                [M, N, dof] = size(motions);
                
                if isempty(bInds)
                    bInds = ones(Np, 1);
                end
                
                if isempty(mInds)
                    mInds = ones(dof, 1);
                end
                
                p = cell(M, N);

                for m = 1:M
                    for n = 1:N
                        pmn = zeros(Np, 1);
                        for o = 1:Np
                            if points(o,3) <= 0
                                for q = 1:dof
                                    if bInds(o) == mInds(q)
                                        del = motions(m, n, q)*funcs(q).Evaluate(points(o,:));
                                        pmn(o) = pmn(o) + -rho*g*del(3);  
                                    end
                                end
                            end
                        end
                        p{m, n} = pmn;
                    end
                end
            else
                p = wf.pressure(type);
            end
            
            pgeo = cell(size(p));
            N = numel(p);
            
            for n = 1:N
                pgeon = PanelGeo(wf.bodyGeo);
                pgeon.Values = p{n};
                pgeo{n} = pgeon;
            end
        end
                
        % Elevation
        function [eta] = Elevation(wf, type)
            error('Elevation not implemented for BodySurfWaveField');
        end
        
        % Velocity
        function [velgeo] = Velocity(wf, type)
            vel = wf.velocity(type);
            velgeo = cell(size(vel));
            N = numel(vel);
            
            for n = 1:N
                velgeon = PanelGeo(wf.bodyGeo);
                velgeon.Values = vel{n};
                velgeo{n} = velgeon;
            end
        end
        
        % SignificantWaveHeight
        function [hs] = SigWaveHeight(wf, type)
            error('SigWaveHeight not implemented for BodySurfWaveField');
        end
        
        % Gets spectrums at the points closest to the desired points
        function [specs, actPoints] = Spectra(wf, varargin)
            error('Spectra not implemented for BodySurfWaveField');
        end
        
        % EnergyFlux
        function [fluxgeo] = EnergyFlux(wf, surf, varargin)
            fluxgeo = wf.bodyGeo;
            
            pnts = fluxgeo.Centroids;
            nrms = fluxgeo.Normals;
            ars = fluxgeo.Areas;
            surf = ControlSurface(pnts, nrms, ars);
            
            flux = wf.energyFlux(surf, varargin{:}, 'PerPoint');

            fluxgeo = cell(size(flux));
            N = numel(flux);
            
            for n = 1:N
                fluxgeon = PanelGeo(wf.bodyGeo);
                fluxgeon.Values = flux{n};
                fluxgeo{n} = fluxgeon;
            end
        end
        
        function [Fe] = WritePressureFile(wf, folder, name, time, iT, A, bodyInds, bodyNames, deltaPoints, compForce)
            % if time less than 1, compute hydrostatic pressure
            if nargin < 4
                time = 0;
            end
            
            if nargin < 5
                iT = 1;
            end
            
            if nargin < 6
                A = 1;
            end
            
            if nargin < 7
                bodyInds = [];
            end
            
            if nargin < 8
                bodyNames = [];
            end
            
            if nargin < 9
                deltaPoints = [];
            end
            if isempty(deltaPoints)
                deltaPoints = [0 0 0];
            end
            
            Fe = [];
            if nargin < 10
                compForce = false;
            end
            
            if isempty(iT)
                iT = 1;
            end
                                    
            geo = wf.bodyGeo;
            cents = geo.Centroids;
            
            T = wf.T(iT);
                        
            pressH0 = wf.Pressure('Hydrostatic');
            pressH0 = pressH0{1}.Values;
            pressHs = wf.Pressure('MotionHydrostatic');
            pressHs = pressHs{iT}.Values;
            pressHd = wf.Pressure('Total');
            pressHd = pressHd{iT}.Values;
            
            Np = geo.Count;
            
            if isempty(bodyInds)
                bodyInds = 1;
            end
            
            if compForce
                funcs = wf.motionFuncs;
                minds = wf.MotionBodyInds;
                Nf = length(funcs);
                Fe = cell(length(time), length(bodyInds));
            end
                        
            for m = 1:length(time)
                
                if time(m) < 0
                    At = 1;
                    
                    tstr = 'pressHydrostatic';
                else
                    phase = wrapTo2Pi(time(m)/T*2*pi);
                    At = A*exp(1i*phase);
                    
                    tstr = num2str(time(m), '%4.2f');
                    tstr(tstr == '.') = '-';
                    tstr = ['t' tstr];
                end
                
                for nb = 1:length(bodyInds)
                    
                    if compForce
                        finds = zeros(Np, 1);
                    end
                    
                    if isempty(bodyNames)
                        bstr = ['bod' num2str(nb)];
                    else
                        bstr = bodyNames{nb};
                    end
                                        
                    if time(m) < 0
                        fileName = [folder '\' name '_' tstr '_' bstr '.csv'];
                    else
                        fileName = [folder '\' name '_' tstr '_press_' bstr '.csv'];
                    end

                    header = {'x (m)', 'y (m)', 'z (m)', 'pressure (Pa)'};

                    N = sum(wf.bodyInds == bodyInds(nb));
                    M = zeros(N, 4);
                    np = 1;
                    for n = 1:Np
                        if wf.bodyInds(n) == bodyInds(nb)
                            M(np, 1:3) = cents(n,:) + deltaPoints;
                            if time(m) < 0
                                M(np, 4) = pressH0(n);
                            else
                                M(np, 4) = pressH0(n) + real(At*(pressHs(n) + pressHd(n)));
                            end

                            np = np+1;
                            if compForce
                                finds(n) = 1;
                            end
                        end
                    end

                    fid = fopen(fileName, 'w');
                    fprintf(fid, '%s,%s,%s,%s\n', header{1}, header{2}, header{3}, header{4});

                    for n = 1:N
                        fprintf(fid, '%7.3f,%7.3f,%7.3f,%8.0f\n', M(n,1), M(n,2), M(n,3), M(n,4));
                    end
                    fclose(fid);
                    
                    if compForce
                        Fe{m, nb} = zeros(3, 1);
                        finds = logical(finds);
                        P = zeros(size(finds));
                        P(finds) = M(:,4);
                        for n = 1:Nf
                            if minds(n) == bodyInds(nb)
                                if bodyInds(nb) == 2
                                    Fe{m, nb}(n-3) = BodySurfWaveField.ComputeForceFromPress(geo, funcs(n), P, finds);
                                else
                                    Fe{m, nb}(n) = BodySurfWaveField.ComputeForceFromPress(geo, funcs(n), P, finds);
                                end
                            end
                        end
                    end
                end
            end
           
        end
        
        function [forceHs0, forceHsM, forceHd, forceHs0e, forceHsMe, forceHde] ...
                = CheckPressureForce(wf, mass, dofSur, dofHea, dofPit, hydroForces, type, iT)
            
            if nargin < 7
                type = 'Total';
            end
            if nargin < 8
                iT = 1;
            end
            
            xi = wf.BodyMotions;
            
            pressHs0 = [];
            pressHsM = [];
            pressHd = [];
            pressTot = [];
            
            if strcmpi(type, 'Hydrostatic') || strcmpi(type, 'Total')
                press = wf.Pressure('Hydrostatic');
                pressHs0 = press{1};
            end
            
            if strcmpi(type, 'MotionHydrostatic') || strcmpi(type, 'Total')
                press = wf.Pressure('MotionHydrostatic');
                pressHsM = press{iT};
            end
            
            if strcmpi(type, 'Dynamic') || strcmpi(type, 'Total')
                press = wf.Pressure('Total');
                pressHd = press{iT};
            end
            
%             if strcmpi(type, 'Total')
%                 pressTot = pressHs0.Values + pressHd.Values + pressHsM.Values;
%             end
                        
            geo = wf.bodyGeo;
            funcs = wf.motionFuncs;
            bInds = wf.bodyInds;
            mInds = wf.motBodyInds;
            Nf = length(funcs);
            
            forceHs0 = zeros(Nf, 1);
            forceHsM = zeros(Nf, 1);
            forceHd = zeros(Nf, 1);
            forceTot = zeros(Nf, 1);
            
            for n = 1:Nf
                
                inds = mInds(n) == bInds;
                
                if ~isempty(pressHs0)
                    forceHs0(n) = BodySurfWaveField.ComputeForceFromPress(geo, funcs(n), pressHs0.Values, inds);
                end
                if ~isempty(pressHsM)
                    forceHsM(n) = BodySurfWaveField.ComputeForceFromPress(geo, funcs(n), pressHsM.Values, inds);
                end
                if ~isempty(pressHd)
                    forceHd(n) = BodySurfWaveField.ComputeForceFromPress(geo, funcs(n), pressHd.Values, inds);    
                end
                if ~isempty(pressTot)
                    forceTot(n) = BodySurfWaveField.ComputeForceFromPress(geo, funcs(n), pressTot, inds);
                end
            end
            
            xi = squeeze(xi(iT, 1, :));
            g = IWaves.G;
            
            if isrow(mass)
                mass = mass';
            end
            
            forceHs0e = zeros(size(forceHs0));
            forceHs0e(dofHea) = -mass*g;
            
            forceHsMe = zeros(size(forceHsM));
            forceHsMe(dofSur) = g*mass.*xi(dofPit);
            forceHsMe = forceHsMe + hydroForces.C*xi;
            
            Fex = squeeze(hydroForces.Fex(iT,1,:));
            A = squeeze(hydroForces.A(iT,:,:));
            B = squeeze(hydroForces.B(iT,:,:));
            w = 2*pi/hydroForces.T(iT);
            forceHde = -(Fex + w^2*A*xi - 1i*w*B*xi);
            
            if strcmpi(type, 'Hydrostatic')
                force = forceHs0;
                forceExp = forceHs0e;
            elseif strcmpi(type, 'MotionHydrostatic')
                force = forceHsM;
                forceExp = forceHsMe;
            elseif strcmpi(type, 'Dynamic')
                force = forceHd;
                forceExp = forceHde;
            elseif strcmpi(type, 'Total');
                force = forceHs0 + forceHsM + forceHd;
                forceExp = forceHs0e + forceHsMe + forceHde;
            end
            
            err = abs(force-forceExp)./abs(forceExp)*100;
            
            fprintf('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
            
            for n = 1:length(force)
                if ~isa(funcs(n), 'ZeroMotionFunc')
                    fprintf('\nBody: %i,\t%s\nForce \t\t\t= %4.2e,\nForce Expected \t= %4.2e\nError \t\t\t= %4.1f\n\n',...
                        mInds(n), class(funcs(n)), abs(force(n)), abs(forceExp(n)), err(n))
                end
            end
        end
    end
    
    methods (Static)
        function [force] = ComputeForceFromPress(geo, func, press, inds)            
            N = geo.Count;
                        
            if nargin < 4
                inds = ones(N, 1);
            else
                N = sum(inds);
            end
            
            cents = geo.Centroids;
            norms = geo.Normals;
            areas = geo.Areas;
            
            cents = cents(inds, :);
            norms = norms(inds, :);
            areas = areas(inds);
            press = press(inds);
            
            normMot = zeros(N, 3);
            dotNorm = zeros(N, 1);
            
            for n = 1:N
                normMot(n, :) = func.Evaluate(cents(n,:));
                dotNorm(n) = dot(normMot(n,:), norms(n,:));
            end

            force = sum(dotNorm.*areas.*press);
        end
    end
end