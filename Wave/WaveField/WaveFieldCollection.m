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
classdef WaveFieldCollection < IWaveField & handle
    % A wave field provides information (pressure, elevation, velocity,
    % significant wave height, spectrums) at discrete points in space
    % specified on either a grid (array) or at arbitrary (x,y,z) points.  
    % 
    % A WaveFieldCollection contains a set of WaveFields that all have the 
    % same periods. The WaveFieldCollection is organized by another 
    % parameter called an index (such as wave direction). These indicies 
    % are given a type (i.e. name), which must be a string and should 
    % describe the collection. The indices themselves do not need to be 
    % numeric. For example, they could be a cell array of strings that 
    % describes each wave field.
    %
    % All the WaveFields in the collection must be equivalent, which means,
    % they have: values at the same spatial locations, the same water 
    % density, the same periods, and the same water depth.
    %
    % Outputs (pressure, elevation, etc.) from the WaveFieldCollection are 
    % organized by wave period then the collection index, because wave 
    % period is the primary orgainzation parameter for wave fields. That 
    % is, because in linear wave theory, waves at different periods do not 
    % interact, while waves at the same period may interact (depending on 
    % how they are defined).
    
    properties (Access = private)
        wavefields;
        inds;
        nInd;
        type;
        collOfColl;
    end

    properties (Dependent)
        Indices;                       % Indicies that organize the wave fields (do not have to be numberical)
        CollType;                  % Name (type) of indices - describes how the wave field collection is oganized
        WaveFields;                 % Array of the internal non-directional wave fields
        WFcount;                    % Number of WaveFields
    end
    
    methods
        
        % Constructor
        function [wf] = WaveFieldCollection(waveFields, varargin)
            % The only necessary argument is an array of WaveFields.
            % Optionally, these two arguments can be supplied together 
            %   'type' - A descriptor (string) of the type of collection
            %   'indices' - A regular (or cell) array of values that
            %   correspond to each WaveField. Must be the same length as
            %   the WaveFields
            
            if (nargin > 0)  
                if (~isempty(varargin))
                    typ = varargin{1};
                    indices = varargin{2};
                    
                    if (~ischar(typ)) 
                        error('WaveFieldCollection name must be a string');
                    end
                    
                    wf.type = typ;
                    wf.inds = indices;
                    
                    wf.nInd = length(indices);
                    
                    if (wf.nInd ~= length(waveFields))
                        error('Number of wave fields is not the wave as the number of indices.');
                    end
                else
                    wf.type = '';
                    wf.nInd = length(waveFields);
                    wf.inds = 1:wf.nInd;
                end                               
                
                wf.hasvel = true;

                wf.wavefields = WaveField.empty(wf.nInd, 0);

                wavefld1 = waveFields(1);

                for n = 1:wf.nInd
                    if (~isa(waveFields(n), 'IWaveField'))
                        error('Input is not a IWaveField object');
                    end
                    if (waveFields(n) ~= wavefld1)
                        error('Not all wave fields are equal');
                    end
                    if (~waveFields(n).hasvel)
                        wf.hasvel = false;
                    end
                end
                
                wf.isarray = wavefld1.IsArray;
                wf.rho = wavefld1.Rho;
                wf.g = wavefld1.G;
                wf.t = wavefld1.T;
                wf.nT = length(wf.t);
                wf.h = wavefld1.H;
                
                if (wavefld1.IsArray)
                    [X Y] = wavefld1.FieldPoints;
                    wf.x = X;
                    wf.y = Y;
                else
                    pts = wavefld1.FieldPoints;
                    wf.points = pts;
                end

                if (isa(waveFields(1), 'WaveFieldCollection'))
                    wf.collOfColl = true;
                else
                    wf.collOfColl = false;
                end
                
                wf.wavefields = waveFields;
            end
        end
                
        function [bet] = get.Indices(wf)
            % Values of the indices that define the WaveFieldCollection
            bet = wf.inds;
        end
        
        function [indNam] = get.CollType(wf)
            % The name of index that organized the collection. e.g.
            % 'Direction'
            indNam = wf.type;
        end
        
        function [wvflds] = get.WaveFields(wf)
            % The WaveFields in the collection
            wvflds = wf.wavefields;
        end
        
        function [wfcnt] = get.WFcount(wf)
            % The number of WaveFields in the collection
            wfcnt = wf.nInd;
        end
        
        % Pressure
        function [pr] = Pressure(wf, varargin)
            % complex pressure values at wave field points.  Optional input
            % of 'Surface' available which returns only pressure values at
            % z = 0
            
            pr = cell(wf.GetTotalSize);
            for n = 1:wf.nInd
                % all these extra :,: is to accomodate collections of
                % collections that are computed recursively. (Technically,
                % only 10 levels of collections are allowed based on the
                % number of :
                pr(:,n,:,:,:,:,:,:,:,:,:,:) = wf.wavefields(n).Pressure(varargin{:});
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
                siz = wf.GetTotalSize;
                v = cell([siz, 3]);

                for n = 1:wf.nInd
                    % all these extra :,: is to accomodate collections of
                % collections that are computed recursively. (Technically,
                % only 10 levels of collections are allowed based on the
                % number of :
                    v(:,n,:,:,:,:,:,:,:,:,:,:,:) = wf.wavefields(n).Velocity(varargin{:});
                end
            else
                v = [];
            end
        end
        
        function [hs] = SigWaveHeight(wf, varargin)
            % spectral significant wave height at surface field points. 
            % Significant wave height is returned for each wave field 
            % separately unless the optional argument, 'Merge',  is used. 
            % 'Merge' combines the significant wave height computation for 
            % all of the wave fields. This would be done, for example, when
            % the WaveFieldCollection is a collection of plane waves at
            % different directions that are considered to represent a given
            % sea.
            
            merge = checkOptions({'Merge'}, varargin);
            
            eta = wf.Elevation;
            
            if (merge)
                if (wf.isarray) 
                    hs = zeros(size(wf.x));
                else
                    hs = zeros(1,size(wf.points, 1));
                end
            else
                hs = cell(wf.nInd, 1);
            end

            for n = 1:wf.nInd
                if (~merge)
                    if (wf.isarray) 
                        hs{n} = zeros(size(wf.x));
                    else
                        hs{n} = zeros(1,size(wf.points, 1));
                    end
                end
                for m = 1:wf.nT
                    thisHs = 0.5*abs(eta{m,n}).^2;
                    if (merge)
                        hs(:,:) = hs(:,:) + thisHs;
                    else
                        hs{n} = hs{n} + thisHs;
                    end
                end
                if (~merge)
                    hs{n} = 4*sqrt(hs{n});
                end
            end
            
            if (merge)
                hs = 4*sqrt(hs);
            end
        end
        
        function [specs, actPoints] = Spectra(wf, varargin)
            % Returns spectra at either all points in the wave field or specified
            % points (optional argument 'Points', points (Nx2). 
            % Spectra are returned for each wave field separately unless 
            % the optional argument, 'Merge',  is used. 
            % 'Merge' produces a directional spectrum for 
            % all of the wave fields. This would be done, for example, when
            % the WaveFieldCollection is a collection of plane waves at
            % different directions that are considered to represent a given
            % sea.
            
            [opts, args] = checkOptions({{'Points', 1}, 'Merge'}, varargin);
            
            usePts = opts(1);
            pts = [];
            if (usePts)
                pts = args{1};
            end
            
            merge = opts(2);
            
            if (merge)
                eta = wf.Elevation;

                f = 1./wf.T;
                df = computeDelta(f);
                db = 1;
                beta = wf.inds;
                if (length(beta) > 1)
                    db = computeDelta(beta);
                    [indb] = find (db < 0);
                    db(indb) = db(indb) + 180;
                end

                [Df, Db] = meshgrid(df, db);

                if (usePts)
                    if (~wf.isarray)
                        error('Point option for Spectrum is only valid in an array wave field.');
                    end

                    [npts, col] = size(pts);
                    if (col ~= 2)
                        error('points must be an Nx2 array');
                    end

                    actPoints = NaN(npts, 2);

                    specs(npts, 1) = WaveSpectrum;

                    for n = 1:npts
                        xi = pts(n,1);
                        yi = pts(n,2);
                        [xo ix yo iy] = wf.FindClosestPoint(xi, yi);
                        actPoints(n,:) = [xo yo];

                        a_ = zeros(wf.nInd, wf.nT);
                        for m = 1:wf.nT
                            for o = 1:wf.nInd
                                a_(o, m) = abs(eta{m, o}(iy, ix));
                            end
                        end

                        S = (0.5*a_.^2./Df./Db).';

                        specs(n) = WaveSpectrum(S, f, beta);
                    end
                else
                    if (wf.isarray)
                        xv = wf.x(1,:);
                        yv = wf.y(:,1);

                        specs(length(yv), length(xv)) = WaveSpectrum;

                        for n = 1:length(yv)
                            for m = 1:length(xv)

                                a_ = zeros(wf.nInd, wf.nT);
                                for o = 1:wf.nT
                                    for q = 1:wf.nInd
                                        a_(q, o) = abs(eta{o, q}(n, m));
                                    end
                                end

                                S = (0.5*a_.^2./Df./Db);

                                specs(n, m) = WaveSpectrum(S, f, beta);
                            end
                        end
                    else
                        npts = length(eta{1,1});
                        specs(npts,1) = WaveSpectrum;

                        for n = 1:npts
                            a_ = zeros(wf.nInd, wf.nT);
                                for o = 1:wf.nT
                                    for q = 1:wf.nInd
                                        a_(q, o) = abs(eta{o, q}(n));
                                    end
                                end

                            S = (0.5*a_.^2./Df./Db);

                            specs(n) = WaveSpectrum(S, f, beta);
                        end
                    end

                    actPoints = [];
                end
            else
                specs = cell(wf.nInd, 1);
                
                for n = 1:wf.nInd
                    if (usePts)
                        specs{n} = wf.wavefields(n).Spectra('Points', pts);
                    else
                        specs{n} = wf.wavefields(n).Spectra;
                    end
                end
            end
        end
        
        % EnergyFlux
        function [flux] = EnergyFlux(wf, surf, varargin)
            if (~wf.hasvel)
                error('Cannot compute flux. Wave field does not contain velocity components.');
            end
            
            flux = cell(wf.nT, wf.nInd);
            for n = 1:wf.nInd
                fluxB = wf.wavefields(n).EnergyFlux(surf, varargin{:});
                for m = 1:wf.nT
                    flux{m,n} = fluxB{m};
                end
            end
        end
        
        function [] = RemoveGeometries(wf, bodies)
            for n = 1:wf.nInd
                wf.wavefields(n).RemoveGeometries(bodies);
            end
        end
        
        function [siz] = GetTotalSize(wf)
            if (wf.collOfColl)
                thatsiz = wf.wavefields(1).GetTotalSize;
                siz = [wf.nT, wf.nInd, thatsiz(2:end)];
            else
                siz = [wf.nT, wf.nInd];
            end
        end
              
        % Overloaded equality operator
        function [areEq] = eq(wfa, wfb)   
            if (~isa(wfa, 'IWaveField') || ~isa(wfb, 'IWaveField'))
                error('Each argument must be a IWaveField');
            end
            
            if (isa(wfb, 'WaveFieldCollection'))
                bColl = true;
            else
                bColl = false;
            end
            
            areEq = 1;
                
            if (bColl)  
                if (wfa.wavefields(1) ~= wfb.wavefields(1))
                    areEq = 0;
                    return;
                end
                    
                if (~strcmp(wfa.type, wfb.type))
                    areEq = 0;
                    return;
                end

                % indicies
                if (wfa.nInd ~= wfb.nInd)
                    areEq = 0;
                    return;
                end

                indCell = false;
                if (iscell(wfa.inds))
                    indCell = true;
                end

                if (iscell(wfb.inds) && ~indCell)
                    areEq = 0;
                    return;
                end

                for n = 1:wfa.nInd
                    if (indCell)
                        if (wfa.inds{n} ~= wfb.inds{n})
                            areEq = 0;
                            return;
                        end
                    else
                        if (wfa.inds(n) ~= wfb.inds(n))
                            areEq = 0;
                            return;
                        end
                    end
                end
            else
                if (wfa.wavefields(1) ~= wfb)
                    areEq = 0;
                    return;
                end
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
               error('Wave fields are not equivalent.  Cannot add.');
           end
           
           if (isa(wfb, 'WaveFieldCollection'))
               bColl = true;
           else
               bColl = false;
           end
           
           if (bColl)               
               for n = 1:wfa.nInd
                   waveflds(n) = wfa.wavefields(n) + wfb.wavefields(n);
               end
               
               wfout = WaveFieldCollection(waveflds, wfa.type, wfa.inds);
           else
               waveflds(wfa.nInd, 1) = WaveField;
               for n = 1:wfa.nInd
                   waveflds(n) = wfa.wavefields(n) + wfb;
               end
               
               wfout = WaveFieldCollection(waveflds, wfa.type, wfa.inds);
            end
        end
        
        % Overloaded unary minus operator
        function [wfout] = uminus(wfin)
            waveflds(wfin.nInd, 1) = WaveField;
            
            for n = 1:wfin.nInd
                waveflds(n) = -wfin.wavefields(n);
            end
            
            wfout = WaveFieldCollection(waveflds, wfin.type, wfin.inds);
        end
        
        % Overloaded unary minus operator
        function [wfout] = minus(wfa, wfb)
            wfb = -wfb;
            wfout = wfa + wfb;
        end
        
        % Overloaded matrix multiplication operator
        function [wfout] = mtimes(a, b)
            if (isnumeric(a) && isa(b, 'WaveFieldCollection'))
                num = a;
                wfin = b;
            elseif (isnumeric(b) && isa(a, 'WaveFieldCollection'))
                num = b;
                wfin = a;
            else
                error('Multiplication only defined between a number and a WaveFieldCollection.')
            end
            
            ni = wfin.nInd;
            wavefdsin = wfin.wavefields;
            
            if (length(num) > 1)
                if any(size(num) ~= wfin.GetTotalSize)
                    error('Numeric argument must either be a scalar or a matrix of same total size as the WaveFieldCollection');
                end

                for n = 1:ni
                    waveflds(n) = squeeze(num(:,n,:,:,:,:,:,:,:,:,:,:))*wavefdsin(1); % index is 1 bc setting the wavefldsin = [] reduces the number of wave fields
                    wavefdsin(1) = [];
                end
            else
                for n = 1:ni
                    waveflds(n) = num*wavefdsin(1);
                    wavefdsin(1) = [];
                end
            end
                 
            wfout = WaveFieldCollection(waveflds, wfin.type, wfin.inds);
        end
        
        % Overloaded multiplication operator
        function [wfout] = times(a, b)
            wfout = a*b;
        end
    end
end