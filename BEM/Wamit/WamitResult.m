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
classdef WamitResult < handle
    % Reads Wamit output files
    %
    % It can be set up by using the WamitRunConditions used to create the 
    % run as the constructor argument.  
    
    % Otherwise, it requires the path to the folder in which the output 
    % files are located and it 
    
    properties (Access = private)
        rho;
        g;
        h;
        folder;
        runName;
        floatingbodies;
        t;
        nT;
        beta;
        nB;
        dof;
        fieldPoints;
        fieldArray;
        solveForces;
        solveField;
        solveBody;
        hydroForces;
        waveBody;
        wavePoints;
        waveArray;
        readVelocity;
        hasBeenRead;
    end

    properties (Dependent)
        Rho;                % Fluid density
        G;                  % Gravitational constant
        H;                  % Depth
        Folder;             % Folder location of run (string)
        RunName;            % Name of Run (string)
        FloatingBodies;     % FloatingBodies (plural because of multiple bodies)
        T;                  % Periods (s)
        Beta;               % Directions (deg)
        HydroForces;        % Hydrodynamic forces computed by WAMIT
        WavePoints;         % Wave field at field points
        WaveArray;          % Wave field at array locations
        WaveBody;           % Wave field on surface of body
    end

    methods

        function [result] = WamitResult(runCondition)
            % Constructor
            if (nargin == 0)
                result.folder = ' ';
                result.runName = [];
            else 
                if (isa(runCondition, 'WamitRunCondition'))
                    result.rho = runCondition.Rho;
                    result.folder = runCondition.Folder;
                    result.runName = runCondition.RunName;
                    result.floatingbodies = runCondition.FloatingBodies;
                    result.dof = HydroBodyComp.GetDoF(result.floatingbodies);
                    result.t = runCondition.T;
                    result.nT = length(result.t);
                    result.beta = runCondition.Beta;
                    result.nB = length(result.beta);
                    result.h = runCondition.H;
                    result.solveForces = runCondition.SolveRad;
                    result.solveBody = runCondition.ComputeBodyPoints;
                    result.fieldPoints = runCondition.FieldPoints;
                    result.fieldArray = runCondition.FieldArray;
                    if (~isempty(result.fieldPoints) || ~isempty(result.fieldArray))
                        result.solveField = true;
                    end
                    result.readVelocity = runCondition.ComputeVelocity;
                else
                    error('Constructor input must be of type WamitRunCondition')
                end
            end
            result.hasBeenRead = 0;
        end
        
        function [fol] = get.Folder(result)
            % Get the folder location where the output files are found
            fol = result.folder;
        end
        function [] = set.Folder(result, fol)
            % Set the folder location where the output files are found
            if (ischar(fol))
                result.folder = fol;
            else
                error('Folder must be a string');
            end
            result.hasBeenRead = 0;
        end
        
        function [rn] = get.RunName(result)
            % Get the name of the run - if Wamit_RunCondition is used,
            % .cfg, .pot, and .frc all have the same file name which is the
            % run name.  Otherwise, the run name is given after the .frc
            % file, which corresponds to the names of the output files.
            rn = result.runName;
        end
        function [] = set.RunName(result, rn)
            % Set the name of the run - if Wamit_RunCondition is used,
            % .cfg, .pot, and .frc all have the same file name which is the
            % run name.  Otherwise, the run name is given after the .frc
            % file, which corresponds to the names of the output files.
            if (ischar(rn))
                result.runName = rn;
            else
                error('RunName must be a string');
            end
            result.hasBeenRead = 0;
        end

        function [fbs] = get.FloatingBodies(result)
            % Get the array of floating bodies
            fbs = result.floatingbodies;
        end
        function [] = set.FloatingBodies(result, fbs)
            % Set the array of floating bodies
            for n = 1:length(fbs)
                if (~isa(fbs(n), 'FloatingBody'))
                    error('All Floating Bodies must be of type FloatingBody');
                end
            end
            
            result.floatingbodies = fbs;
            result.hasBeenRead = 0;
        end
        
        function [rh] = get.Rho(result)
            % The water density
            if (result.hasBeenRead)
                rh = result.rho;
            else
                rh = [];
            end
        end
        
        function [g_] = get.G(result)
            % The gravitational constant
            if (result.hasBeenRead)
                g_ = result.g;
            else
                g_ = [];
            end
        end
        
        function [t_] = get.T(result)
            % The wave periods
            if (result.hasBeenRead)
                t_ = result.t;
            else
                t_ = [];
            end
        end

        function [bet] = get.Beta(result)
            % The wave headings
            if (result.hasBeenRead)
                bet = result.beta;
            else
                bet = [];
            end
        end
               
        function [h_] = get.H(result)
            % The water depth
            if (result.hasBeenRead)
                h_ = result.h;
            else
                h_ = [];
            end
        end
        
        function [hf] = get.HydroForces(result)
            % The HydroForces object
            if (result.hasBeenRead)
                hf = result.hydroForces;
            else
                hf = [];
            end
        end
        
        function [wp] = get.WavePoints(result)
            % Wave field points, not in an array
            if (result.hasBeenRead)
                wp = result.wavePoints;
            else
                wp = [];
            end
        end
        
        function [wa] = get.WaveArray(result)
            % Wave field points in an array
            if (result.hasBeenRead)
                wa = result.waveArray;
            else
                wa = [];
            end
        end
        
        function [wb] = get.WaveBody(result)
            % Wave field points in an array
            if (result.hasBeenRead)
                wb = result.waveBody;
            else
                wb = [];
            end
        end

        function [] = ReadResult(result, varargin)
            % Reads the Wamit results

            useSing = false;
            removeBodies = false;
            for n = 1:length(varargin)
                if (strcmp(varargin{n}, 'UseSingle'))
                    useSing = true;
                end
                if (strcmp(varargin{n}, 'RemoveBodies'))
                    removeBodies = true;
                end
            end
                        
            % If no run name is given, get it from .frc in fnames.wam
            if (isempty(result.runName))
                fid = fopen([result.folder '\fnames.wam']);
                if (fid == -1)
                    error('Run name could not be located');
                end
                name = [];
                while(~feof(fid))
                    line = fgetl(fid);
                    [pa, name, ext] = fileparts(line);
                    if (strcmp(ext, '.frc'))
                        result.runName = name;
                        break;
                    end
                end
                if (isempty(name))
                    error('Run name could not be located');
                else
                    result.runName = name;
                end
                fclose(fid);
            end
            
            % Gravity, depth and density
            fid = fopen([result.folder '\' result.runName '.out']);
            if (fid == -1)
                error('No WAMIT results for this run!');
            end
            while(~feof(fid))
                line = fgetl(fid);
                line = strtrim(line);
                if (length(line) > 7)
                    if (strcmp(line(1:8), 'Gravity:'))
                        g_ = str2num(strtrim(line(9:20)));
                        result.g = g_;
                        
                        line = fgetl(fid);                        
                        
                        h_ = str2num(strtrim(line(20:29)));
                        if (isempty(h_))
                            h_ = Inf;
                        end
                        if (isempty(result.h))
                            result.h = h_;
                        else
                            if (isinf(h_))
                                if (~isinf(result.h))
                                    error('Depth read in is not the same as the expected depth');
                                end
                            else
                                if (abs((h_ - result.h)/h_) > 2e-3)
                                    error('Depth read in is not the same as the expected depth');
                                end
                            end
                        end
                        
                        rh = str2num(strtrim(line(56:65)));
                        if (isempty(result.rho))
                            result.rho = rh;
                        else
                            if (abs((rh - result.rho)/rh) > 2e-3)
                                error('Density read in is not the same as the expected density');
                            end
                        end
                        
                        break;
                    end
                    
                end
            end
            fclose(fid);
                        
            fullpath = result.folder;
            
            if (isempty(result.solveForces))
                result.solveForces = false;
                fid1 = fopen([result.folder '\' result.runName '.1']);
                fid2 = fopen([result.folder '\' result.runName '.2']);
                
                if ((fid1 ~= -1) && (fid2 ~= -1))
                    result.solveForces = true;
                end

                if (fid1 ~= -1)
                    fclose(fid1);
                end
                
                if (fid2 ~= -1)
                    fclose(fid2);
                end
            end
            
            if (result.solveForces)                
                % Added Mass and Damping 
                [a_, b_, t_, modes] = Wamit_read1(fullpath, result.runName, result.rho);
                result.dof = length(modes);
                
                % Check periods
                if (isempty(result.t))
                    result.t = t_;
                    result.nT = length(result.t);
                else
                    if (length(t_) ~= result.nT)
                        error('Periods read in are not the same as the expected periods');
                    end

                    for n = 1:result.nT
                        if (abs((t_(n) - result.t(n))/t_(n)) > 2e-3)
                            error('Periods read in are not the same as the expected periods');
                        end
                    end
                end
                
                % Hydrostatic matrix
                c1 = Wamit_readHst(fullpath, result.runName, result.rho, result.g);
                c_ = zeros(result.dof, result.dof);
                for m = 1:result.dof
                    for n = 1:result.dof
                        c_(m,n) = c1(modes(m), modes(n));
                    end
                end

                % Excitation force
                [f, t_, bet] = Wamit_read2(fullpath, result.runName, result.rho, result.g);
                
                % Check directions
                if (isempty(result.beta))
                    result.beta = bet;
                    result.nB = length(result.beta);
                else
                    if (length(bet) ~= result.nB)
                        error('Directions read in are not the same as the expected directions')
                    end

                    % Check directions
    %                 for n = 1:result.nB
    %                     if (abs(bet(n) - result.beta(n)) > 1e-4)
    %                         error('Directions read in are not the same as the expected directions')
    %                     end
    %                 end
                end

                result.hydroForces = HydroForces(result.t, result.beta, a_, b_, c_, f, result.h, result.rho);
            end
            
            if (isempty(result.solveBody))
                result.solveBody = false;
                fid1 = fopen([result.folder '\' result.runName '.5p']);
                fid2 = fopen([result.folder '\' result.runName '.5vx']);
                
                if (fid1 ~= -1)
                    result.solveBody = true;
                    fclose(fid1);
                end
                
                if (fid2 ~= -1)
                    result.readVelocity = true;
                    fclose(fid2);
                end
                
            end
            
            if (result.solveBody)
                
                if (isempty(result.floatingbodies))
                    error('WamitResult must contain a floating body to compute body points');
                end
                
                if (length(result.floatingbodies) > 1)
                    error('WamitResult can only contain a single floating body to computer body points.');
                end
                bodyGeo = result.floatingbodies(1).PanelGeo;
                % bodyGeo is given in body coordinates, where some of to
                % make sure it is submerged
                bodyGeo.Translate([0, 0, result.floatingbodies(1).Zpos]);

                [P_rad, P_diff, bodCenters] = Wamit_readNum5p(fullpath, result.runName, result.floatingbodies(1).GeoFile, length(result.t), length(result.beta), result.dof, result.rho, result.g, useSing);
                
                if (result.readVelocity)
                    [V_rad, V_diff] = Wamit_readNum5v(fullpath, result.runName, result.floatingbodies(1).GeoFile, result.t, length(result.beta), result.dof, result.g);
                end
                
                nBodPoints = size(bodCenters, 1);
                if (nBodPoints ~= bodyGeo.Count)
                    error('The number of body points must be the same as the number of panels.');
                end
                
                centsBG = bodyGeo.Centroids;
                
                for n = 1:nBodPoints
                    if any(abs(centsBG(n,:) - bodCenters(n,:)) > 1e-3*[1 1 1])
                        %error('The body points and wave field points must be the same');
                    end
                end
                
                p_rad_pts = P_rad(:, :, 1:nBodPoints);
                p_diff_pts = P_diff(:, :, 1:nBodPoints);
                
                    if (result.readVelocity)
                        v_rad_pts = V_rad(:, :, :, 1:nBodPoints);
                        v_diff_pts = V_diff(:, :, :, 1:nBodPoints);
                    else
                        v_rad_pts = [];
                        v_diff_pts = [];
                    end
                    
                    for n = 1:result.nB
                        thisP = zeros(result.nT, nBodPoints);
                        if (useSing)
                            thisP = single(thisP);
                        end
                        if (result.readVelocity)
                            thisV = zeros(result.nT, 3, nBodPoints);
                        else
                            thisV = [];
                        end

                        for m = 1:result.nT
                            thisP(m, :) = squeeze(p_diff_pts(m, n, :));
                            if (result.readVelocity)
                                thisV(m, :, :) = squeeze(v_diff_pts(m, n, :, :));
                            end
                        end

                        dwfn = WaveField(result.rho, result.g, result.h, result.t, thisP, thisV, 0, bodCenters);
                        wcs = PlaneWaves(ones(size(result.t)), result.t, result.beta(n)*ones(size(result.t)), result.h);
                        iwfn = PlaneWaveField(result.rho, wcs, false, bodCenters);
                        iwfs(n) = iwfn;
                        swfs(n) = dwfn - iwfn;
                    end
                    
                    iwavefield = WaveFieldCollection(iwfs, 'Direction', result.beta); 
                    swavefield = WaveFieldCollection(swfs, 'Direction', result.beta);

                    rwfs(result.dof, 1) = WaveField;
                    for n = 1:result.dof
                        thisP = zeros(result.nT, nBodPoints);
                        if (useSing)
                            thisP = single(thisP);
                        end
                        if (result.readVelocity)
                            thisV = zeros(result.nT, 3, nBodPoints);
                        else
                            thisV = [];
                        end

                        for m = 1:result.nT
                            thisP(m, :) = squeeze(p_rad_pts(m, n, :));
                            if (result.readVelocity)
                                thisV(m, :, :) = squeeze(v_rad_pts(m, n, :, :));
                            end
                        end

                        rwfs(n) = WaveField(result.rho, result.g, result.h, result.t, thisP, thisV, 0, bodCenters);
                    end
                    rwavefield = WaveFieldCollection(rwfs, 'MotionIndex', (1:result.dof));

                    result.waveBody = BodySurfWaveField(bodyGeo, iwavefield, swavefield, rwavefield);
            end
            
            if (isempty(result.solveField))
                result.solveField = false;
                fid1 = fopen([result.folder '\' result.runName '.6p']);
                fid2 = fopen([result.folder '\' result.runName '.6vx']);
                
                if (fid1 ~= -1)
                    result.solveField = true;
                    fclose(fid1);
                end
                
                if (fid2 ~= -1)
                    result.readVelocity = true;
                    fclose(fid2);
                end
                
            end

            % New read wave field points code for 7.0
            % Read in field points and arrays
            if (result.solveField)
                [P_rad, P_diff, points] = Wamit_readNum6p(fullpath, result.runName, length(result.t), length(result.beta), result.dof, result.rho, result.g, useSing);
                
                if (result.readVelocity)
                    [V_rad, V_diff] = Wamit_readNum6v(fullpath, result.runName, result.t, length(result.beta), result.dof, result.g);
                end
                
                if (~isempty(result.fieldPoints))
                    nPoints = size(result.fieldPoints, 1);
                elseif(~isempty(result.fieldArray))
                    nPoints = 0;
                else
                    nPoints = size(P_rad, 3);
                    data = importdata([result.folder '\' result.runName '.fpt']);
                    data = data.data;
                    result.fieldPoints = squeeze(data(:,2:4));
                end
                
                % Field Points
                if (nPoints > 0)
                    p_rad_pts = P_rad(:, :, 1:nPoints);
                    p_diff_pts = P_diff(:, :, 1:nPoints);

                    if (result.readVelocity)
                        v_rad_pts = V_rad(:, :, :, 1:nPoints);
                        v_diff_pts = V_diff(:, :, :, 1:nPoints);
                    else
                        v_rad_pts = [];
                        v_diff_pts = [];
                    end
                    
                    for n = 1:result.nB
                        thisP = zeros(result.nT, nPoints);
                        if (useSing)
                            thisP = single(thisP);
                        end
                        if (result.readVelocity)
                            thisV = zeros(result.nT, 3, nPoints);
                        else
                            thisV = [];
                        end

                        for m = 1:result.nT
                            thisP(m, :) = squeeze(p_diff_pts(m, n, :));
                            if (result.readVelocity)
                                thisV(m, :, :) = squeeze(v_diff_pts(m, n, :, :));
                            end
                        end

                        dwfn = WaveField(result.rho, result.g, result.h, result.t, thisP, thisV, 0, result.fieldPoints);
                        wcs = PlaneWaves(ones(size(result.t)), result.t, result.beta(n)*ones(size(result.t)), result.h);
                        iwfn = PlaneWaveField(result.rho, wcs, false, result.fieldPoints);
                        iwfs(n) = iwfn;
                        swfs(n) = dwfn - iwfn;
                    end
                    
                    iwavefield = WaveFieldCollection(iwfs, 'Direction', result.beta); 
                    swavefield = WaveFieldCollection(swfs, 'Direction', result.beta);

                    rwfs(result.dof, 1) = WaveField;
                    for n = 1:result.dof
                        thisP = zeros(result.nT, nPoints);
                        if (useSing)
                            thisP = single(thisP);
                        end
                        if (result.readVelocity)
                            thisV = zeros(result.nT, 3, nPoints);
                        else
                            thisV = [];
                        end

                        for m = 1:result.nT
                            thisP(m, :) = squeeze(p_rad_pts(m, n, :));
                            if (result.readVelocity)
                                thisV(m, :, :) = squeeze(v_rad_pts(m, n, :, :));
                            end
                        end

                        rwfs(n) = WaveField(result.rho, result.g, result.h, result.t, thisP, thisV, 0, result.fieldPoints);
                    end
                    rwavefield = WaveFieldCollection(rwfs, 'MotionIndex', (1:result.dof));

                    result.wavePoints = FBWaveField(iwavefield, swavefield, rwavefield);
                end

                % Field Array
                if (~isempty(result.fieldArray))
                    points = points((nPoints + 1):end, :);
                    if (nPoints ~= 0)
                        P_rad = P_rad(:, :, (nPoints + 1):end);
                        P_diff = P_diff(:, :, (nPoints + 1):end);
                    end

                    [P_rad, P_diff] = WamitResult.makeRectangular(P_rad, P_diff, points, useSing);

                    if (result.readVelocity)
                        if (nPoints ~= 0)
                            V_rad = V_rad(:, :, :, (nPoints + 1):end);
                            V_diff = V_diff(:, :, :, (nPoints + 1):end);
                        end

                        [V_rad, V_diff] = WamitResult.makeRectangular(V_rad, V_diff, points);
                    else
                        V_rad = [];
                        V_diff = [];
                    end

                    [X, Y] = result.fieldArray.GetArrayPoints;

                    removeBodies = true;

                    if (removeBodies)
                        val = 0.05*result.rho*result.g;
                        indBodies = (abs(squeeze(P_diff(1, 1, :, :))) < val);
                    end

                    for n = 1:result.nB
                        thisP = zeros([result.nT, size(X)]);
                        if (useSing)
                            thisP = single(thisP);
                        end
                        if (result.readVelocity)
                            thisV = zeros([result.nT, 3, size(X)]);
                        else
                            thisV = [];
                        end

                        for m = 1:result.nT
                            thisP(m, :, :) = squeeze(P_diff(m, n, :, :));
                            if (result.readVelocity)
                                thisV(m, :, :, :) = squeeze(V_diff(m, n, :, :, :));
                            end

                            if (removeBodies)
                                thisP(m, indBodies) = NaN;
                                if (result.readVelocity)
                                    for o = 1:3
                                        thisV(m, o, indBodies) = NaN;
                                    end
                                end
                            end
                        end

                        dwfn = WaveField(result.rho, result.g, result.h, result.t, thisP, thisV, 1, X, Y);
                        wcs = PlaneWaves(ones(size(result.t)), result.t, result.beta(n)*ones(size(result.t)), result.h);
                        iwfn = PlaneWaveField(result.rho, wcs, 1, X, Y);
                        iwfs(n) = iwfn;
                        swfs(n) = dwfn - iwfn;
                    end

                    iwavefield = WaveFieldCollection(iwfs, 'Direction', result.beta); 
                    swavefield = WaveFieldCollection(swfs, 'Direction', result.beta);

                    rwfs(result.dof, 1) = WaveField;

                    for n = 1:result.dof
                        thisP = zeros([result.nT, size(X)]);
                        if (useSing)
                            thisP = single(thisP);
                        end
                        if (result.readVelocity)
                            thisV = zeros([result.nT, 3, size(X)]);
                        else
                            thisV = [];
                        end

                        for m = 1:result.nT
                            thisP(m, :, :) = squeeze(P_rad(m, n, :, :));
                            if (result.readVelocity)
                                thisV(m, :, :, :) = squeeze(V_rad(m, n, :, :, :));
                            end
                        end

                        if (removeBodies)
                            thisP(m, indBodies) = NaN;
                            if (result.readVelocity)
                                for o = 1:3
                                    thisV(m, o, indBodies) = NaN;
                                end
                            end
                        end

                        rwfs(n) = WaveField(result.rho, result.g, result.h, result.t, thisP, thisV, 1, X, Y);
                    end
                    
                    rwavefield = WaveFieldCollection(rwfs, 'MotionIndex', (1:result.dof));

                    result.waveArray = FBWaveField(iwavefield, swavefield, rwavefield);
                end
            end
            
            result.hasBeenRead = 1;
        end
    end
    
    % Static Methods    
    methods (Static, Access = private)
        %{
        % makeRectangular method for 6.0
        function [wf2, X, Y] = makeRectangular(wf1, points, varargin)
            useSing = 0;
            
            if (~isempty(varargin))
                useSing = varargin{1};
            end

            x = unique(points(:,1));
            y = unique(points(:,2));
            [X, Y] = meshgrid(x,y);
            Nx = length(x);
            Ny = length(y);

            s = size(wf1);

            if length(s) == 3
                [Nper Ndof buffer] = size(wf1);

                if (useSing)
                    wf2 = single(zeros(Nper, Ndof, Ny, Nx));
                else
                    wf2 = zeros(Nper, Ndof, Ny, Nx);
                end

                for m = 1:Ny
                    for n = 1:Nx
                        [I, F] = mode([find(points(:,1) == x(n)); find(points(:,2) == y(m))]);
                        if F==2
                            wf2(:, :, m, n) = wf1(:, :, I);
                        end
                    end
                end
            elseif length(s) == 4
                [Nper Ndof buffer buffer] = size(wf1);

                wf2 = zeros(Nper, Ndof, 3, Ny, Nx);

                for m = 1:Ny
                    for n = 1:Nx
                        [I, F] = mode([find(points(:,1) == x(n)); find(points(:,2) == y(m))]);
                        if F==2
                            wf2(:, :, :, m, n) = wf1(:, :, :, I);
                        end
                    end
                end
            end
        end
        %}
        
        % makeRectangular function for 7.0
        function [rad2, diff2, X, Y] = makeRectangular(rad, diff, points, varargin)
            useSing = 0;
            
            if (~isempty(varargin))
                useSing = varargin{1};
            end

            x = unique(points(:,1));
            y = unique(points(:,2));
            [X, Y] = meshgrid(x,y);
            Nx = length(x);
            Ny = length(y);

            s = size(rad);

            if length(s) == 3
                [Nper Ndof buffer] = size(rad);
                [buffer Nbeta buffer] = size(diff);

                if (useSing)
                    rad2 = single(zeros(Nper, Ndof, Ny, Nx));
                    diff2 = single(zeros(Nper, Nbeta, Ny, Nx));
                else
                    rad2 = zeros(Nper, Ndof, Ny, Nx);
                    diff2 = zeros(Nper, Nbeta, Ny, Nx);
                end

                for m = 1:Ny
                    for n = 1:Nx
                        [I, F] = mode([find(points(:,1) == x(n)); find(points(:,2) == y(m))]);
                        if F==2
                            rad2(:, :, m, n) = rad(:, :, I);
                            diff2(:, :, m, n) = diff(:, :, I);
                        end
                    end
                end
            elseif length(s) == 4
                [Nper Ndof buffer buffer] = size(rad);
                [buffer Nbeta buffer buffer] = size(diff);

                rad2 = zeros(Nper, Ndof, 3, Ny, Nx);
                diff2 = zeros(Nper, Nbeta, 3, Ny, Nx);

                for m = 1:Ny
                    for n = 1:Nx
                        [I, F] = mode([find(points(:,1) == x(n)); find(points(:,2) == y(m))]);
                        if F==2
                            rad2(:, :, :, m, n) = rad(:, :, :, I);
                            diff2(:, :, :, m, n) = diff(:, :, :, I);
                        end
                    end
                end
            end
        end
    end
end