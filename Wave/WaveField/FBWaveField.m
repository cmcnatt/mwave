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
classdef FBWaveField < IWaveField & handle
    % A wave field provides information (pressure, elevation, velocity,
    % significant wave height, spectrums) at discrete points in space
    % specified on either a grid (array) or at arbitrary (x,y,z) points.  
    %
    % The FBWaveField (FB = 'FloatingBody') is a WaveField that handles the
    % floating body situation. It contains:
    % 
    %   - Incident wave field(s) - must be IIncWaveField. One or more. 
    %       Multiple wave fields must be supplied as an IncidentWaveField, 
    %       which is a WaveFieldCollection that implements IIncWaveField. 
    %   - Scattered wave field(s) - can be any type of WaveField. One or 
    %       more. Multiple wave fields must be supplied as a 
    %       WaveFieldCollection. The number of scattered wave fields must 
    %       equal the number of incident wave fields. 
    %   - Radiatted wave field(s) - can be any type of WaveField. 
    %       Zero or more. The number of radiated wave fields must equal the
    %       total number of degrees of freedom of the floating body or 
    %       floating body array. Zero radiated wave fields means the body 
    %       or bodies are held fixed. Supplied at object initialization.
    %   - Body motion(s) - an array of complex amplitudes which are the 
    %       amplitude and phase of the body motions in each degree of 
    %       freedom of the body or array due to each incident wave. The 
    %       size of the body motins array must be 
    %       [Nperiod, Nincident, Ndof]. Optional argument at object 
    %       initialization, and can be supplied or changed through the 
    %       'BodyMotions' property after object initialization. 

    
    properties (Access = private)
        iwavefield;
        swavefield;
        rwavefield;
        isPlaneInc;
        indices;
        nInc;
        dof;
        motions;
        a;
    end

    properties (Dependent)
        IncWaveVals;
        BodyMotions;                % Complex amplitudes of the body motions
        DoF;                        % Degrees of freedom
        IncidentWaveAmps;           % Wave amplitude of incoming wave Nper x Nbeta matrix of complex numbers (amplitude and phase)
    end
    
    methods
        
        % Constructor
        function [wf] = FBWaveField(iwave, swave, varargin)
                                    
            if (isa(iwave, 'WaveFieldCollection'))
                if (~isa(swave, 'WaveFieldCollection'))
                    error('If the incident wave is a WaveFieldCollection, the scattered wave must be a WaveFieldCollection also.');
                end
                
                if (iwave ~= swave)
                    error('The incident and scattered waves must be equivalent');
                end     
                
                wf.nInc = swave.WFcount;
                wf.indices = swave.Indices;
                swf1 = swave.WaveFields(1);
            else
                if (isa(swave, 'WaveFieldCollection'))
                    error('If the incident wave is not a WaveFieldCollection, then the scattered wave cannot be WaveFieldCollection either.');
                end
                if (length(iwave) ~= 1)
                    error('Mulitple incident waves must be supplied as a WaveFieldCollection');
                end
                
                if (length(swave) ~= 1)
                    error('Mulitple scattered waves must be supplied as a WaveFieldCollection');
                end
                
                wf.nInc = 1;
                wf.indices = [];
                swf1 = swave;
            end
            
            wf.hasvel = true;
            
            if (~iwave.hasvel)
                wf.hasvel = false;
            end
                        
            if (~swave.hasvel)
                wf.hasvel = false;
            end
            
            if (~isempty(varargin))
                rwave = varargin{1};
                if (length(rwave) > 1)
                    error('Mulitple radiated waves must be supplied as a WaveFieldCollection');
                else
                    if (isa(rwave, 'WaveFieldCollection'))
                        wf.dof = rwave.WFcount;
                    else
                        if (~isa(rwave, 'WaveField'))
                            error('Radiated wave field is not a WaveField object');
                        end
                        wf.dof = 1;
                    end
                end
                
                if (~rwave.hasvel)
                    wf.hasvel = false;
                end
                
                if (wf.dof > 1)
                    rwf1 = rwave.WaveFields(1);
                else
                    rwf1 = rwave;
                end
                if (swf1 ~= rwf1)
                    error('Not all wave fields are equal');
                end
                
                wf.rwavefield = rwave;
            else
                wf.dof = 0;    
                wf.rwavefield = [];
            end
            
            [opts, args] = checkOptions({{'BodyMotions', 1}}, varargin);
            
            if (opts(1))
                bodyMotions = args{1};
                [bm1, bm2, bm3] = size(bodyMotions);
                if (bm3 ~= wf.dof)
                    error('The number of body motions must be the same as the number of radiated wave fields.');
                end
                            
                wf.motions = bodyMotions;
            end
            
            wf.isarray = swave.IsArray;
            wf.rho = swave.Rho;
            wf.g = swave.G;
            wf.t = swave.T;
            wf.nT = length(wf.t);
            wf.h = swave.H;
                                   
            if (wf.isarray)
                [X, Y] = swave.FieldPoints;
                wf.x = X;
                wf.y = Y;
            else
                pts = swave.FieldPoints;
                wf.points = pts;
            end
            
            wf.iwavefield = iwave;
            wf.swavefield = swave;
        end
        
        % Beta
        function [bet] = get.IncWaveVals(wf)
            bet = wf.indices;
        end
        
        % BodyMotions
        function [mot] = get.BodyMotions(wf)
            mot = wf.motions;
        end
        function [] = set.BodyMotions(wf, mot)
            if (~isempty(mot))
                [nt, nb, df] = size(mot);
                if ((nt ~= wf.nT) || (nb ~= wf.nInc) || (df ~= wf.dof))
                    error('The size of the body motions array must be nT x nInc x DoF');
                end
            end
            wf.motions = mot;
        end
        
        % DoF
        function [df] = get.DoF(wf)
            df = wf.dof;
        end
        
        % IncidentWaves
        function [a_] = get.IncidentWaveAmps(wf)
            a_ = wf.a;
        end
        function [] = set.IncidentWaveAmps(wf, a_)
            [row, col] = size(a_);
            
            if (row ~= wf.nT && col ~= wf.nInc)
                error('The incident wave array must be a matrix of size (number of periods) by (number of incident waves)');
            end
            
            wf.a = a_;
        end
                
        % Pressure
        function [p] = Pressure(wf, type)
            p = wf.pressure(type);
        end
                
        % Elevation
        function [eta] = Elevation(wf, type)
            p = wf.pressure(type);

            eta = cellfun(@(x) x./(wf.rho*wf.g), p, 'UniformOutput', false);
        end
        
        % Velocity
        function [vel] = Velocity(wf, type)
            vel = wf.velocity(type);
        end
        
        % SignificantWaveHeight
        function [hs] = SigWaveHeight(wf, type)
            if (nargin < 1)
                type = 'Total';
            end
            
            thiswf = wf.getWF(type);
            hs = thiswf.SigWaveHeight('Merge');
%             if (strcmp(wf.iwavefield.CollType, 'Direction'))
%                 hs = thiswf.SigWaveHeight('Merge');
%             else
%                 hs = thiswf.SigWaveHeight;
%             end
        end
        
        % Gets spectrums at the points closest to the desired points
        function [specs, actPoints] = Spectra(wf, varargin)
            type = 'Total';
            usePts = ' ';
            pts = [];
                        
            for n = 1:length(varargin)
                if (strcmp(varargin{n}, 'Points'))
                    usePts = 'Points';
                    pts = varargin{n+1};
                elseif (strcmp(varargin{n}, 'Radiated'))
                    type = 'Radiated';
                elseif (strcmp(varargin{n}, 'Diffracted'))
                    type = 'Diffracted';
                elseif (strcmp(varargin{n}, 'Incident'))
                    type = 'Incident';
                elseif (strcmp(varargin{n}, 'Scattered'))
                    type = 'Scattered';
                end
            end
            
            thiswf = wf.getWF(type);
            if (strcmp(wf.iwavefield.CollType, 'Direction'))
                [specs, actPoints] = thiswf.Spectra('Merge', usePts, pts);
            else
                [specs, actPoints] = thiswf.Spectra(usePts, pts);
            end
        end
        
        % EnergyFlux
        function [flux] = EnergyFlux(wf, surf, varargin)
            flux = wf.energyFlux(surf, varargin{:});
        end
        
        function [] = RemoveGeometries(wf, bodies)
            wf.swavefield.RemoveGeometries(bodies);
            wf.rwavefield.RemoveGeometries(bodies);
        end
              
        % Overloaded equality operator
        function [areEq] = eq(wfa, wfb)
            if (~isa(wfa, 'FBWaveField') || ~isa(wfb, 'FBWaveField'))
                error('Each argument must be a FBWaveField');
            end
            
            areEq = 1;
            
            % dof
            if any(wfa.dof ~= wfb.dof)
                areEq = 0;
                return;
            end
            
            % iwavefield - checks the points, periods and directions, ets
            if (wfa.iwavefield ~= wfb.iwaveflied)
                areEq = 0;
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
            error('Addition operator not defined for FBWaveField');
        end
        
        % Overloaded unary minus operator
        function [wfout] = uminus(wfin)
            error('Negation operator not defined for FBWaveField');
        end
        
        % Overloaded unary minus operator
        function [wfout] = minus(wfa, wfb)
            error('Subtraction operator not defined for FBWaveField');
        end
        
        % Overloaded matrix multiplication operator
        function [wfout] = mtimes(a, b)
            error('Multiplication operator not defined for FBWaveField');
        end
        
        % Overloaded multiplication operator
        function [wfout] = times(a, b)
            wfout = a*b;
        end
    end
    
     methods (Access = protected)
         
        function [p] = pressure(wf, type)
            if (nargin < 1)
                type = 'Total';
            end
            
            thiswf = wf.getWF(type);
            p = thiswf.Pressure;
        end
        
        function [vel] = velocity(wf, type)
            if (wf.hasvel)            
                if (nargin < 1)
                    type = 'Total';
                end
                
                thiswf = wf.getWF(type);
                vel = thiswf.Velocity;
            else
                vel = [];
            end
        end
     
        function [flux] = energyFlux(wf, surf, varargin)
            if (~wf.hasvel)
                error('Cannot compute flux. Wave field does not contain velocity components.');
            end
                        
            opt = [];
            if (isempty(varargin))
                type = 'Total';
            else
                type = varargin{1};
                if (length(varargin) > 1)
                    opt = varargin{2};
                end
            end
            
            thiswf = wf.getWF(type);
            flux  = thiswf.EnergyFlux(surf, opt);
        end
         
         function [thiswf] = getWF(wf, type)
              switch (type)
                case 'Incident'
                    if (isempty(wf.a))
                        thiswf = wf.iwavefield;
                    else
                        thiswf = wf.a*wf.iwavefield;
                    end
                case 'Radiated'
                    if (isempty(wf.motions))
                        thiswf = wf.rwavefield;
                    else
                        rwfsMot(wf.nInc, 1) = WaveFieldCollection;
                        
                        for m = 1:wf.nInc
                            xi = squeeze(wf.motions(:,m,:));
                            if (wf.nT == 1)
                                xi = xi.';
                            end
                            rwfsMot(m) = xi.*wf.rwavefield;
                        end
                        
                        thiswf = WaveFieldCollection(rwfsMot, wf.iwavefield.CollType, wf.iwavefield.Indices);
                    end
                case 'Diffracted'
                    if (isempty(wf.a))
                        thiswf = wf.iwavefield + wf.swavefield;
                    else
                        thiswf = wf.a*(wf.iwavefield + wf.swavefield);
                    end
                case 'Scattered'
                    if (isempty(wf.a))
                        thiswf = wf.swavefield;
                    else
                        thiswf = wf.a*wf.swavefield;
                    end
                case 'Total'                    
                    if (wf.dof > 0)
                        if isempty(wf.motions)
                            error('Cannot return Total wave field. Body motions are empty.');
                        end
                        
                        rwfs = wf.rwavefield.WaveFields;
                        for m = 1:wf.nInc
                            rwfsmn = ZeroWaveField(rwfs(1));
                            for n = 1:wf.dof
                                rwfsmn = rwfsmn + squeeze(wf.motions(:,m,n)).*rwfs(n);
                            end
                            rwfsMot(m) = rwfsmn;
                        end
                        
                        radwf = WaveFieldCollection(rwfsMot, wf.iwavefield.CollType, wf.iwavefield.Indices);
                        if (isempty(wf.a))
                            thiswf = wf.iwavefield + wf.swavefield + radwf;
                        else
                            thiswf = wf.a*(wf.iwavefield + wf.swavefield) + radwf;
                        end
                        
                    else
                        if isempty(wf.a)
                            thiswf = wf.iwavefield + wf.swavefield;
                        else
                            thiswf = wf.a*(wf.iwavefield + wf.swavefield);
                        end
                    end
                otherwise
                    error('Wave type must be ''Incident'', ''Radiated'', ''Diffracted'', ''Scattered'', or ''Total''');
              end
         end
    end
end