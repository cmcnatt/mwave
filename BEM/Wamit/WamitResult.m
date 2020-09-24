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
classdef WamitResult < IBemResult
    % Reads Wamit output files
    %
    % It can be set up by using the WamitRunConditions used to create the 
    % run as the constructor argument.  
    
    % Otherwise, it requires the path to the folder in which the output 
    % files are located and it 
    
    properties (Access = private)
        runName;
        solveRad;
        solveDiff;
        solveBody;
        waveBody;
        readVelocity;
        compFK;
        fieldPointsAsBody;
        bodyPansWithPress;
        bodyPanInds;
        compDrift;
        driftOption;
        iReadRAO;
        spikeFreqs;
    end

    properties (Dependent)
        RunName;            % Name of Run (string)
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
                    result.dof = FreqDomComp.GetDoF(result.floatingbodies);
                    result.t = runCondition.T;
                    result.nT = length(result.t);
                    result.beta = runCondition.Beta;
                    result.nB = length(result.beta);
                    result.h = runCondition.H;
                    result.solveRad = runCondition.SolveRad;
                    result.solveDiff = runCondition.SolveDiff;
                    result.solveBody = runCondition.ComputeBodyPoints;
                    result.fieldPoints = runCondition.FieldPoints;
                    result.fieldPointsAsBody = ~isempty(runCondition.FieldPointAsBodyDistance);
                    result.bodyPansWithPress = runCondition.BodyPanelsWithPress;
                    result.bodyPanInds = runCondition.BodyPanelIndices;
                    result.fieldArray = runCondition.FieldArray;
                    result.cylArray = runCondition.CylArray;
                    result.compFK = runCondition.ComputeFK;
                    result.compDrift = runCondition.CompDrift;
                    result.driftOption = runCondition.DriftOption;
                    result.iReadRAO = runCondition.IReadRAO;
                    result.spikeFreqs = runCondition.SpikeFreqs;
                    if (~isempty(result.fieldPoints) || ~isempty(result.fieldArray) || ~isempty(result.cylArray))
                        result.solveField = true;
                    end
                    result.readVelocity = runCondition.ComputeVelocity;
                else
                    error('Constructor input must be of type WamitRunCondition')
                end
            end
            result.errLog = [];
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

            % read error log
            fid = fopen([result.folder '\errorp.log']);
            nerr = 0;
            while(~feof(fid))
                line = fgetl(fid);
                if (~isempty(line))
                    if (~isempty(line(3:4)))
                        if (strcmp(line(3:4), 'No'))
                            line = fgetl(fid);
                            line = fgetl(fid);
                            nerr = nerr + 1;
                            result.errLog{nerr} = strtrim(line);
                        end
                    end
                end
            end
            
            fclose(fid);
            
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
                warning('No WAMIT results for this run!');
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
            
            if (isempty(result.solveRad) || isempty(result.solveDiff))
                result.solveRad = false;
                result.solveDiff = false;
                fid1 = fopen([result.folder '\' result.runName '.1']);
                fid2 = fopen([result.folder '\' result.runName '.2']);
                fid3 = fopen([result.folder '\' result.runName '.3']);
                fid4 = fopen([result.folder '\' result.runName '.2fk']);
                fid5 = fopen([result.folder '\' result.runName '.3fk']);
                
                if fid1 ~= -1
                    result.solveRad = true;
                end
                
                if (fid2 ~= -1) && (fid3 ~= -1) && (fid4 ~= -1) && (fid5 ~= -1)
                    result.solveDiff = true;
                end

                if (fid1 ~= -1)
                    fclose(fid1);
                end
                
                if (fid2 ~= -1)
                    fclose(fid2);
                end
                
                if (fid3 ~= -1)
                    fclose(fid3);
                end
                
                if (fid4 ~= -1)
                    fclose(fid4);
                    result.compFK = true;
                end
                
                if (fid4 ~= -1)
                    fclose(fid4);
                    result.compFK = true;
                end
            end
            
            % Defualt Values
            a_ = [];
            b_ = [];
            c_ = [];
            ainf = [];
            a0 = [];
            f = [];
            f_fk = [];
            
            if (result.solveRad)                
                % Added Mass and Damping 
                if ~isempty(result.spikeFreqs{1,1}) || ~isempty(result.spikeFreqs{2,1}) % i.e. if spike interpolation frequencies have been provided
                    [a_, b_, t_, modes, ainf, a0] = Wamit_read1(fullpath, result.runName, result.rho, 'removeSpikes', result.spikeFreqs);
                else
                    [a_, b_, t_, modes, ainf, a0] = Wamit_read1(fullpath, result.runName, result.rho);
                end
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

            end
            
            if (result.solveDiff)
                % Excitation force
                if result.compFK
                    if ~isempty(result.spikeFreqs{1,1}) || ~isempty(result.spikeFreqs{2,1}) % i.e. if spike interpolation frequencies have been provided
                        [f_fk] = Wamit_read23(fullpath, result.runName, result.rho, result.g, 'fk', 'removeSpikes', result.spikeFreqs);
                        [f_sc, t_, bet] = Wamit_read23(fullpath, result.runName, result.rho, result.g, 'sc', 'removeSpikes', result.spikeFreqs);
                    else
                        [f_fk] = Wamit_read23(fullpath, result.runName, result.rho, result.g, 'fk');
                        [f_sc, t_, bet] = Wamit_read23(fullpath, result.runName, result.rho, result.g, 'sc');
                    end
                    
                    f = f_fk + f_sc;
                else
                    f_fk = [];
                    if ~isempty(result.spikeFreqs{1,1}) || ~isempty(result.spikeFreqs{2,1}) % i.e. if spike interpolation frequencies have been provided
                        [f, t_, bet] = Wamit_read23(fullpath, result.runName, result.rho, result.g, 'removeSpikes', result.spikeFreqs);
                    else
                        [f, t_, bet] = Wamit_read23(fullpath, result.runName, result.rho, result.g);
                    end
                end
                
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
            end
            
            if result.solveRad || result.solveDiff
                result.hydroForces = FreqDomForces(result.t, result.beta, a_, b_, c_, f, result.h, result.rho, f_fk, a0, ainf);
            end
            
            if (result.iReadRAO ~= 1) % i.e. do the following unless only options 1-4 have been computed.
            if (isempty(result.solveBody))
                result.solveBody = false;
                if result.fieldPointsAsBody
                    fid1 = fopen([result.folder '\' result.runName '.6p']);
                    fid2 = fopen([result.folder '\' result.runName '.6vx']);
                else
                    fid1 = fopen([result.folder '\' result.runName '.5p']);
                    fid2 = fopen([result.folder '\' result.runName '.5vx']);
                end
                
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
                
                if result.fieldPointsAsBody 
                    result.solveField = false;
                end
                
                opts = {};
                if ~isempty(result.floatingbodies(1).CompGeoFile)
                    opts{length(opts)+1} = 'bpo';
                end
                
                if result.floatingbodies(1).WamILowHi
                    opts{length(opts)+1} = 'ihi';
                end
                
                cgs = cell(length(result.floatingbodies), 1);
                for n = 1:length(result.floatingbodies)
                    cgs{n} = result.floatingbodies(n).Cg;
                end
                
                if result.fieldPointsAsBody
                    thisDof = result.dof;
                    if ~result.solveRad
                        thisDof = 0;
                    end
                
                    [P_rad, P_diff, centers0] = Wamit_readNum6p(fullpath, result.runName, length(result.t), length(result.beta), thisDof, result.rho, result.g, useSing);
                
                    centers = cell(length(result.floatingbodies), 1);
                    for n = 1:length(result.floatingbodies)
                        thisPan = result.bodyPanInds == n;
                        if sum(thisPan) > 0
                            centers{n} = centers0(thisPan, :);
                        end
                    end
                    
                    if (result.readVelocity)
                        [V_rad, V_diff] = Wamit_readNum6v(fullpath, result.runName, result.t, length(result.beta), thisDof, result.g);
                    end
                else
                    [P_rad, P_diff, centers] = Wamit_readNum5p(fullpath, result.runName, {result.floatingbodies.GeoFile}, length(result.t), length(result.beta), result.dof, result.rho, result.g, cgs, useSing, opts{:});

                    if (result.readVelocity)
                        [V_rad, V_diff] = Wamit_readNum5v(fullpath, result.runName, {result.floatingbodies.GeoFile}, result.t, length(result.beta), result.dof, result.g, useSing, opts);
                    end
                end
                
                ind = 0;
                dofFunc = 0;
                for l = 1:length(result.floatingbodies)

                    motFuncl = result.floatingbodies(l).Modes.MotionFuncs;
                    NmotFunc = length(motFuncl);
                        
                    if ~isempty(centers{l})
                        if result.floatingbodies(l).WamILowHi
                            % high order GDF
                            % For high order, need a CompGeoFile - i.e. a
                            % low order GDF
                            bodyGeo = Wamit_readGdf(fullpath, result.floatingbodies(l).CompGeoFile);
                        else
                            % low order GDF
                            if ~isempty(result.floatingbodies(l).PanelGeo)
                                % if has a PanelGeo us that
                                bodyGeo = result.floatingbodies(l).PanelGeo;
                                % get the wet subset of panels
                                bodyGeo = bodyGeo.GetSubsetPanels(true, false);
                            else
                                bodyGeo = Wamit_readGdf(fullpath, result.floatingbodies(l).GeoFile);
                                % GDF files that are read in are centered
                                % at their center of rotation
                                bodyGeo.Translate(result.floatingbodies(l).CenterRot);
                                % get the wet subset of panels - b/c read
                                % in, interior panels not indicated, so
                                % can't use .GetSubsetPanels, so use
                                % just get the panels where z < 0
                                inds = bodyGeo.Centroids(:,3) < 0;
                                bodyGeo = PanelGeo(bodyGeo.Panels(inds));
                            end
                        end
                        

%                         if isempty(bodyGeo)
%                             if isempty(result.floatingbodies(l).CompGeoFile)
%                                 if ~result.floatingbodies(l).WamILowHi
%                                     meshName = result.floatingbodies(l).GeoFile;
%                                 end
%                             else
%                                 meshName = result.floatingbodies(l).CompGeoFile;
%                             end
%                             bodyGeo = Wamit_readGdf(fullpath, meshName);
%                         else
%                             bodyGeo.Translate([0, 0, result.floatingbodies(1).Zpos]);
%                         end
                        % TODO (line below may need uncommented for higher order meshes)
                        %bodyGeo.Translate(result.floatingbodies(l).Cg);
                        
                        centsBG = bodyGeo.Centroids;
                        [nBodPoints, ~] = size(centsBG);

%                         nBodPoints = size(centers{l}, 1);
%                         if (nBodPoints ~= bodyGeo.Count)
%                             error('The number of body points must be the same as the number of panels.');
%                         end
%                         for n = 1:nBodPoints
%                             if any(abs(centsBG(n,:) - centers{l}(n,:)) > 1e-3*[1 1 1])
%                                 error('The body points and wave field points must be the same');
%                             end
%                         end
%                         bodCenters((ind+1):(ind+nBodPoints),:) = centers{l};

                        panels((ind+1):(ind+nBodPoints)) = bodyGeo.Panels;
                        bodCenters((ind+1):(ind+nBodPoints),:) = centsBG;
                        bodyInds((ind+1):(ind+nBodPoints)) = l;

                        ind = ind + nBodPoints;

                        motionFuncs(dofFunc+1:dofFunc+NmotFunc) = motFuncl;
                    else
                        for q = dofFunc+1:dofFunc+NmotFunc
                            motionFuncs(q) = ZeroMotionFunc;
                        end
                    end
                    
                    motBodyInds(dofFunc+1:dofFunc+NmotFunc) = l;
                    dofFunc = dofFunc + NmotFunc;
                end
                
                bodyGeo = PanelGeo(panels);
                nBodPoints = ind;
                
                if isempty(result.bodyPansWithPress)
                    p_rad_pts = P_rad(:, :, 1:nBodPoints);
                    p_diff_pts = P_diff(:, :, 1:nBodPoints);
                else
                    [M, N, ~] = size(P_rad);
                    p_rad_pts = zeros(M, N, nBodPoints);
                    p_rad_pts(:, :, result.bodyPansWithPress) = P_rad(:, :, :);
                    
                    [M, N, ~] = size(P_diff);
                    p_diff_pts = zeros(M, N, nBodPoints);
                    p_diff_pts(:, :, result.bodyPansWithPress) = P_diff(:, :, :);
                end

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
                result.waveBody = BodySurfWaveField(bodyGeo, iwavefield, swavefield, rwavefield);   % Creates BodySurfWaveField with geometry (PanelGeo), and wave fields (i.e. pressure and velocities)
                result.waveBody.BodyInds = bodyInds;                        % which panels are from which body (if more than one body)
                result.waveBody.MotionFuncs = motionFuncs;                  % Return unit vector for a given motion
                result.waveBody.MotionBodyInds = motBodyInds;
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
                thisDof = result.dof;
                if ~result.solveRad
                    thisDof = 0;
                end                
                
                [P_rad, P_diff, points] = Wamit_readNum6p(fullpath, result.runName, length(result.t), length(result.beta), thisDof, result.rho, result.g, useSing);
                
                if (result.readVelocity)
                    [V_rad, V_diff] = Wamit_readNum6v(fullpath, result.runName, result.t, length(result.beta), thisDof, result.g);
                end
                
                if (~isempty(result.fieldPoints))
                    nPoints = size(result.fieldPoints, 1);
                    fpoints = result.fieldPoints;
                elseif (~isempty(result.cylArray))
                    nPoints = result.cylArray.Ntheta*result.cylArray.Nz;
                    fpoints = result.cylArray.GetPoints(result.h);
                elseif(~isempty(result.fieldArray))
                    nPoints = 0;
                else
                    nPoints = size(P_rad, 3);
                    data = importdata([result.folder '\' result.runName '.fpt']);
                    data = data.data;
                    result.fieldPoints = squeeze(data(:,2:4));
                    fpoints = result.fieldPoints;
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

                        dwfn = WaveField(result.rho, result.g, result.h, result.t, thisP, thisV, 0, fpoints);
                        wcs = PlaneWaves(ones(size(result.t)), result.t, result.beta(n)*ones(size(result.t)), result.h);
                        iwfn = PlaneWaveField(result.rho, wcs, false, fpoints);
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

                        rwfs(n) = WaveField(result.rho, result.g, result.h, result.t, thisP, thisV, 0, fpoints);
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

                    if ~isempty(P_rad)
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
                    else
                        result.waveArray = FBWaveField(iwavefield, swavefield);
                    end
                end
            end
            
            % Read drift forces
            if (result.compDrift)
                result.driftForces = ReadDriftForces(result);
            end
            
            end % if (result.iReadRAO ~= 1)
            
            result.hasBeenRead = 1;
        end
    end
    
    methods (Access = private)
        
        function [df] = ReadDriftForces(result)
            % READDRIFTFORCES Extracts the drift forces from the apropriate 
            % WAMIT output file. Outputs a structure containing all the
            % results. This could be expanded upon later to include further
            % quantities if desired.
            
            driftMethodNames = {'Control surface','Momentum conservation','Pressure integration'};
            
            for n = 1:length(result.driftOption)
                ext = num2str(result.driftOption(n)+6); % Find file extension based on which drift methods have been used
                try
                    data = importdata([result.Folder '\' result.RunName '.' ext]);
                catch
                    error('Unable to extract drift forces from results files.')
                end
                driftData = data.data;
                
                % First, find out how many periods there are and make a vector of them.
                T = unique(driftData(:,1),'stable'); % Want to keep same order so the indices match those of hydro.w
                
                % Second, find out how many directions there are and make a vector for each of them.
                Beta1 = unique(driftData(:,2));
                Beta2 = unique(driftData(:,3));
                
                % Find out which modes are present (including the negative numbers
                % that represent the components of the moment about the moving origin)
                modeIndices = unique(driftData(:,4));
                modeIndices_posOnly = unique(abs(modeIndices)); % The total moments are given by the positively numbered indices
                
                % Dimensionalised mean drift forces for unit amplitude
                % Currently only set up to read drift forces for single
                % frequencies - i.e. the diagonal elements of the QTF
                % matrix.
                mdf = zeros(length(T),max(modeIndices_posOnly),length(Beta1),length(Beta2)); % Initialise mean drift force array
                for j = 1:length(Beta1)
                    for i = 1:length(modeIndices_posOnly) % This depends on the number of bodies present
                        mdf(:,modeIndices_posOnly(i),j,j) = ...
                            (driftData(driftData(:,4)==modeIndices_posOnly(i) & driftData(:,2)==Beta1(j) & driftData(:,3)==Beta2(j),7) + ...
                            1i*driftData(driftData(:,4)==modeIndices_posOnly(i) & driftData(:,2)==Beta1(j) & driftData(:,3)==Beta2(j),8))*9.81*1000;
                    end
                end
                
                df(n).T = T; % Same wave periods have to be used for each method in a given run
                df(n).meanDriftForces = mdf;
                df(n).methodName = driftMethodNames{n};

            end % for n = 1:length(result.driftOption)
            
      %%% NOTE with relation to control surface and pressure integration results:
            %%% This code assumes that we don't care about the division of
            %%% the moments between the static and moving frames, hence
            %%% the data with negative mode indices is ignored here.
            
        end % ReadDriftForces
        
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