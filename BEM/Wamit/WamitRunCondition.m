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
classdef WamitRunCondition < IBemRunCondition
    % Creates input files to a Wamit run
    
    properties (Access = private)
        runName;
        solveDiff;
        solveRad;
        t;
        beta;
        fieldPoints;
        computeVelocity;
        computeBody;
        fieldPointAsBodyDist;
        bodyPansWithPress;
        bodyPanInds;
        scratchPath;
        useridPath;
        ncpu;
        ramGB;
        maxItt;
        useDirect;
        blockItCnt;
        useHaskind;
        computeFK;
        compDrift;
        autoCSF;
        driftOption;
        boxCSF;
        iReadRAO;
        motionRAOs;
    end

    properties (Dependent)
        RunName;            % Name of Run (string) (Nemoh?)
        SolveDiff;          % Solve the diffraction problem
        SolveRad;           % Solve the radiation problem
        T;                  % Periods (s)
        Beta;               % Directions (radians)
        FloatingBodies;     % FloatingBodies (plural because of multiple bodies)
        FieldPoints;        % Arbitrary field points (Nx3 array)
        ComputeBodyPoints;  % Indicates whether pressure and velocity on the body surface will be computed
        FieldPointAsBodyDistance;   % Use the field points .6 instead of body point .5 to compute body surface pressure. 
        BodyPanelsWithPress;        % Indicates whether pressure is computed on the given body panel
        BodyPanelIndices;   % Indices of the bodies that the field points refer to.
        ComputeVelocity;    % Indicates whether velocity should be computed at field points
        WamitPath;          % Path location of wamit.exe
        ScratchPath;        % Path location of wamit scratch folder
        UseridPath;         % Path location of userid.wam license file
        NCPU;               % The number of cpus to be used in the Wamit computation
        RAMGBmax;           % The max RAM to be used in the Wamit computation
        MaxItt;
        UseDirectSolver;
        BlockIterativeSolverCount;
        UseHaskind;
        ComputeFK;
        CompDrift;          % Whether or not to compute drift forces (using control surface(s))
        AutoCSF;            % Whether or not to use automatic control surface creation
        DriftOption;        % Which methods to use to compute the mean drift forces
        BoxCSF;             % Whether to use a 'box' shaped control surface or a cylindrical one
        IReadRAO;           % Whether to use .4 file or user-input RAOs for evaluating options 5-9
        MotionRAOs;         % RAOs computed externally to WAMIT, for computing drift forces with ireadrao option
    end

    methods
        
        function [run] = WamitRunCondition(folder, runName)
            run.rho = 1000;
            run.computeVelocity = false;
            run.computeBody = false;
            run.fieldPointAsBodyDist = [];
            run.bodyPansWithPress = [];
            run.bodyPanInds = [];
            run.solveDiff = true;
            run.solveRad = true;
            run.exePath = 'C:\wamitv7';
            run.scratchPath = '\wamitv7\scratch';
            run.useridPath = '\wamitv7';
            run.ncpu = 1;
            run.ramGB = 2;
            run.maxItt = 35;
            run.useDirect = false;
            run.blockItCnt = [];
            run.useHaskind = false;
            run.computeFK = false;
            run.compDrift = false; % by default do not compute drift forces
            run.driftOption = [];
            run.autoCSF = [];
            run.boxCSF = [];
            run.iReadRAO = 0; % by default, use .4 file to get RAOs for options 5-9
            run.geoFiles = [];
            if (nargin == 0)
                run.folder = ' ';
                run.runName = 'newRun';
            else 
                run.folder = folder;
                run.runName = runName;
            end
        end
        
        function [rn] = get.RunName(run)
            % Get the run name, which is the name that will be given to all
            % of the input files except if they have thier own name, the
            % geometry files
            rn = run.runName;
        end
        function [] = set.RunName(run, val)
            % Set the run name, which is the name that will be given to all
            % of the input files except if they have thier own name, the
            % geometry files
            if (ischar(val))
                run.runName = val;
            else
                error('RunName must be a string');
            end
        end
        
        function [sd] = get.SolveDiff(run)
            % Solve the diffraction problem
            sd = run.solveDiff;
        end
        function [] = set.SolveDiff(run, val)
            % Solve the diffraction problem
            if ((val ~= 1) && (val ~= 0))
                error('value must be boolean, 1 or 0');
            end
            
            run.solveDiff = val;
        end
        
        function [sr] = get.SolveRad(run)
            % Solve the radiation problem
            sr = run.solveRad;
        end
        function [] = set.SolveRad(run, val)
            % Solve the radiation problem
            if ((val ~= 1) && (val ~= 0))
                error('value must be boolean, 1 or 0');
            end
            
            run.solveRad = val;
        end
                
        function [t_] = get.T(run)
            % Get the wave periods to be run in wamit
            t_ = run.t;
        end
        function [] = set.T(run, val)
            % Set the wave periods to be run in wamit
            if (isnumeric(val))
                run.t = val;
            else
                error('The periods must be numbers');
            end
        end
        
        function [bet] = get.Beta(run)
            % Get the wave headings to be run in wamit (radians)
            bet = run.beta;
        end
        function [] = set.Beta(run, val)
            % Set the wave headings to be run in wamit (radians)
            if (isnumeric(val))
                for n = 1:length(val)
                    if(val(n) < -pi || val(n) > 2*pi)
                        error('All wave headings must be between -pi and 2pi radians');
                    end
                end
                run.beta = val;
            else
                error('The wave headings must be numbers');
            end
        end
        
        function [fbs] = get.FloatingBodies(run)
            % Get the array of floating bodies
            fbs = run.floatBods;
        end
        function [] = set.FloatingBodies(run, val)
            % Set the array of floating bodies
            pass = 1;
            genmds = 0;
%             if ((length(val) > 1) && run.computeBody)
%                 error('Cannot compute body points for an array of floating bodies');
%             end
            for n = 1:length(val)
                if (~isa(val(n), 'FloatingBody'))
                    pass = 0;
                end
                if (val(n).Ngen > 0)
                    if (genmds ~= 0)
                        if (genmds ~= val(n).WamIGenMds)
                            error('The Wamit IGENMDS value must be the same for all floating bodes in the array');
                        end
                    else
                        genmds = val(n).WamIGenMds;
                    end
                end
            end
            if (pass)
                run.floatBods = val;
            else
                error('All Floating Bodies must be of type FloatingBody');
            end    
        end
              
        function [fp] = get.FieldPoints(run)
            % Get the field points to be evaluated in wamit
            fp = run.fieldPoints;
        end
        function [] = set.FieldPoints(run, val)
            % Set the field points to be evaluated in wamit
            fail = 1;
            if (ndims(val) == 2)
                [ro co] = size(val);
                if (co == 3)
                    fail = 0;
                end
            end
            if (fail)
                error('Field points must be an Nx3 array');
            end
            run.fieldPoints = val;
        end
                
        function [cb] = get.ComputeBodyPoints(run)
            % Get whether or not pressure and velocity is evaluated on the
            % body
            cb = run.computeBody;
        end
        function [] = set.ComputeBodyPoints(run, val)
            % Set whether or not pressure and velocity is evaluated on the
            % body
            if ((val ~= 1) && (val ~= 0))
                error('value must be boolean, 1 or 0');
            end
            
%             if (length(run.floatBods) > 1)
%                 error('Cannot compute body point for more than one floating body');
%             end
            
            run.computeBody = val;
        end
        
        function [val] = get.FieldPointAsBodyDistance(run)
            % Use the field points .6 instead of body point .5 to compute
            % body surface pressure. The distance is the offset from the
            % body
            val = run.fieldPointAsBodyDist;
        end
        function [] = set.FieldPointAsBodyDistance(run, val)
             % Use the field points .6 instead of body point .5 to compute
            % body surface pressure. The distance is the offset from the
            % body
            
            run.fieldPointAsBodyDist = val;
        end
        
        function [val] = get.BodyPanelsWithPress(run)
            % Indicates on which body panels pressure is computed
            val = run.bodyPansWithPress;
        end
        
        function [val] = get.BodyPanelIndices(run)
            % Indices of the bodies that the field points refer to
            val = run.bodyPanInds;
        end
        
        function [cv] = get.ComputeVelocity(run)
            % Get whether or not velocity is evaluated at field points
            cv = run.computeVelocity;
        end
        function [] = set.ComputeVelocity(run, val)
            % Set whether or not velocity is evaluated at field points
            if ((val ~= 1) && (val ~= 0))
                error('value must be boolean, 1 or 0');
            end
            
            run.computeVelocity = val;
        end
        
        function [sp] = get.ScratchPath(run)
            % Path location of wamit scratch folder
            sp = run.scratchPath;
        end
        function [] = set.ScratchPath(run, val)
            % Path location of wamit scratch folder
            if (ischar(val))
                run.scratchPath = val;
            else
                error('ScratchPath must be a string');
            end
        end
        
        function [idp] = get.UseridPath(run)
            % Path location of userid.wam license file
            idp = run.useridPath;
        end
        function [] = set.UseridPath(run, val)
            % Path location of userid.wam license file
            if (ischar(val))
                run.useridPath = val;
            else
                error('UseridPath must be a string');
            end
        end
        
        function [ncp] = get.NCPU(run)
            % The number of cpus to be used in the Wamit computation
            ncp = run.ncpu;
        end
        function [] = set.NCPU(run, val)
            % The number of CPUs to be used in the Wamit computation
            if (~isInt(val))
                error('The number of CPUs must be an integer');
            end
            if (val < 1)
                error('The number of CPUs must be greater than or equal to 1')
            end
            
            run.ncpu = val;
        end
        
        function [ram] = get.RAMGBmax(run)
            % The max RAM to be used in the Wamit computation
            ram = run.ramGB;
        end
        function [] = set.RAMGBmax(run, val)
            % The max RAM to be used in the Wamit computation
            if (val < 0)
                error('RAMGBmax must be greater than 0');
            end
            
            run.ramGB = val;
        end
        
        function [mi] = get.MaxItt(run)
            % The max number of iterations - default is 35
            mi = run.maxItt;
        end
        function [] = set.MaxItt(run, val)
             % The max number of iterations - default is 35
            if (~isInt(val) || (val <= 0))                
                error('The MaxItt must be an integer greater than 0');
            end
            
            run.maxItt = val;
        end
        
        function [ud] = get.UseDirectSolver(run)
            % Use the direct solver
            ud = run.useDirect;
        end
        function [] = set.UseDirectSolver(run, val)
             % Use the direct solver
            if (~isBool(val))                
                error('UseDirectSolver must be a boolean');
            end
            
            run.useDirect = val;
        end
        
        function [val] = get.BlockIterativeSolverCount(run)
            % Set the block iterative solver count
            val = run.blockItCnt;
        end
        function [] = set.BlockIterativeSolverCount(run, val)
             % Get the block iterative solver count
            if ~isInt(val) || val < 1              
                error('BlockIterativeSolverCount must be a positive integer');
            end
            
            run.blockItCnt = val;
        end
        
        function [val] = get.UseHaskind(run)
            % Use the Haskind relation to compute forces
            val = run.useHaskind;
        end
        function [] = set.UseHaskind(run, val)
             % Use the Haskind relation to compute forces
            if (~isBool(val))                
                error('UseHaskind must be a boolean');
            end
            
            run.useHaskind = val;
        end
        
        function [val] = get.ComputeFK(run)
            % Compute the Froude-Krylov force along with the excitation
            % force
            val = run.computeFK;
        end
        function [] = set.ComputeFK(run, val)
             % Compute the Froude-Krylov force along with the excitation
            % force
            if (~isBool(val))                
                error('ComputeFK must be a boolean');
            end
            
            run.computeFK = val;
        end
        
        function [val] = get.CompDrift(run)
            % Compute the mean drift forces acting on the body/bodies
            val = run.compDrift;
        end
        function [] = set.CompDrift(run, val)
            % Compute the mean drift forces acting on the body/bodies
            if (~isBool(val))                
                error('CompDrift must be a boolean');
            end
            
            run.compDrift = val;
        end
        
        function [val] = get.AutoCSF(run)
            % Compute the mean drift forces acting on the body/bodies
            val = run.autoCSF;
        end
        function [] = set.AutoCSF(run, val)
            % Compute the mean drift forces acting on the body/bodies
            if (~isBool(val))                
                error('AutoCSF must be a boolean');
            end
            
            run.autoCSF = val;
        end
        
        function [val] = get.DriftOption(run)
            % Compute the mean drift forces acting on the body/bodies
            val = run.driftOption;
        end
        function [] = set.DriftOption(run, val)
            % Compute the mean drift forces acting on the body/bodies
            if (~isInt(val))                
                error('DriftOption must be a vector or scalar of only integer values.');
            elseif sum(~ismember(val,[1 2 3]))~=0
                error('Drift options can only take 1, 2 or 3 as possible values.')
            end
            
            run.driftOption = val;
        end
        
        function [val] = get.BoxCSF(run)
            % Use box control surface or a cylindrical one
            val = run.boxCSF;
        end
        function [] = set.BoxCSF(run, val)
            % Use box control surface or a cylindrical one
            if (~isBool(val))                
                error('BoxCSF must be a boolean');
            end
            
            run.boxCSF = val;
        end
        
        function [val] = get.MotionRAOs(run)
            % Pass the motion RAOs in so can be written to .RAO file
            val = run.motionRAOs;
        end
        function [] = set.MotionRAOs(run, val)
            % Pass the motion RAOs in so can be written to .RAO file
            if ~(ndims(val)==3)
                error('MotionRAOs must be a three-dimensional array.')
            end
            % Find total no. of DoFs across all bodies
            nDoFs = 0;
            for i = 1:length(run.floatBods)
                nDoFs = nDoFs + run.floatBods(i).Modes.DoF;
            end
            if (size(val,1)~=length(run.T))                
                error('MotionRAOs must have number of rows equal to the number of wave periods.');
            elseif (size(val,2)~=length(run.Beta))
                error('MotionRAOs must have number of columns equal to the number of wave angles.');
            elseif (size(val,3)~=nDoFs)
                error('MotionRAOs must have same number of columns in direction 3 as the total number of DoFs.')
            end
            
            run.motionRAOs = val; % When the public version, MotionRAOs is set by the user,
                                  % it is essentially relabelled as the
                                  % private version, motionRAOs.
        end
        
        
        function [val] = get.IReadRAO(run)
            % Where to get RAOs from for evaluating options 5-9
            val = run.iReadRAO;
        end
        function [] = set.IReadRAO(run, val)
            % Where to get RAOs from for evaluating options 5-9
            if (~isInt(val))                
                error('IReadRAO must be an integer value.');
            elseif sum(~ismember(val,[0 1 2 3]))~=0
                error('IReadRAO can only take 0, 1, 2 or 3 as possible values.')
            end
            
            run.iReadRAO = val;
        end
        
        
        function [] = CleanRunFolder(run, varargin)
            opts = checkOptions({'ExceptGdf'}, varargin);
            if ~opts(1)
                delete([run.folder '\*.gdf']);
                delete([run.folder '\*.spl']);
                delete([run.folder '\*.ms2']);
            end
            delete([run.folder '\*.1']); 
            delete([run.folder '\*.2']);
            delete([run.folder '\*.3']);
            delete([run.folder '\*.4']);
            delete([run.folder '\*.5*']);
            delete([run.folder '\*.6*']);
            delete([run.folder '\*.7']);
            delete([run.folder '\*.8']);
            delete([run.folder '\*.9']);
            delete([run.folder '\*.frc']);
            delete([run.folder '\*.pot']);
            delete([run.folder '\*.cfg']);
            delete([run.folder '\*.wam']);
            delete([run.folder '\*.bpi']);
            delete([run.folder '\*.rao']);
            delete([run.folder '\*.p2f']);
            delete([run.folder '\*.log']);
            delete([run.folder '\*.txt']);
            delete([run.folder '\*.out']);
            delete([run.folder '\*.num']);
            delete([run.folder '\*.dat']);
            delete([run.folder '\*.pnl']);
            delete([run.folder '\*.hst']);
            delete([run.folder '\*.fpt']);
            delete([run.folder '\*.idf']);
            delete([run.folder '\*.bpo']);
            delete([run.folder '\*.dat']);
            delete([run.folder '\*.csf']);
            delete([run.folder '\*.bat']);
            delete([run.folder '\*.mmx']);
        end
        
        function [] = WriteRun(run, varargin)
            % Writes the Input files (except for the .gdf) for a WAMIT run
            
            [opts] = checkOptions({'NoBatch'}, varargin);
            noBat = opts(1);
            
            if (exist(run.folder, 'file') ~= 7)
                error('Designated run folder does not exist');
            end
                        
            % .gdf and .bpi
            nbody = length(run.floatBods);
            run.geoFiles = cell(1, nbody);
            Np = 0;
            for n = 1:nbody
                if isempty(run.floatBods(n).GeoFile)
                    geoFile = [run.runName num2str(n)];
                    run.writeGdf(run.floatBods(n), geoFile, run.floatBods(n).ISurfPan);
                    run.geoFiles{n} = geoFile;
                    run.floatBods(n).GeoFile = geoFile;
                else
                    run.geoFiles{n} = run.floatBods(n).GeoFile;
                end
                
                if run.computeBody
                    if ~isempty(run.floatBods(n).CompGeoFile)
                        geo = Wamit_readGdf(run.folder, run.floatBods(n).CompGeoFile);
                        geo.Translate(run.floatBods(n).Cg);

                        if isempty(run.fieldPointAsBodyDist)
                            Wamit_writeBpi(run.folder, run.geoFiles{n}, geo.Centroids);
                        else
                            dist = run.fieldPointAsBodyDist;
                            if length(dist) > 1
                                dist = dist(n);
                            end
                            Npn = geo.Count;
                            points(Np+1:Np+Npn,:) = geo.OffsetCentroids(dist);
                            panInds(Np+1:Np+Npn) = n;
                            Np = Np + Npn;
                        end
                    end
                end
            end
            
            if Np > 0 
                if ~isempty(run.fieldPoints)
                    error('FieldPoints not empty. Cannot use field points to compute body surface points.');
                else
                    run.bodyPansWithPress = points(:,3) <= 0;
                    run.bodyPanInds = panInds(run.bodyPansWithPress);
                    run.fieldPoints = points(run.bodyPansWithPress,:);
                end
            end
                        
            if (~noBat)
                filename = [run.folder '\wam_run.bat'];
                fileID = fopen(filename, 'wt');
                
                fprintf(fileID, ['set "fN=' run.runName '"\n\n']);
                fprintf(fileID, 'del %%fN%%.out %%fN%%.1 %%fN%%.2 %%fN%%.3 %%fN%%.4 ');
                fprintf(fileID, '%%fN%%.5 %%fN%%.6p %%fN%%.6vx %%fN%%.6vy %%fN%%.6vz ');
                fprintf(fileID, '%%fN%%.fpt ');
                if run.iReadRAO == 0 %i.e. if running WAMIT for just options 5-9, no need to rerun POTEN module,
                    % so need to retain the p2f file.
                fprintf(fileID, '%%fN%%.p2f ');
                fprintf(fileID, 'errorp.log ');
                end
                fprintf(fileID, 'errorf.log ');
                fprintf(fileID, 'rgkinit.txt rgklog.txt wamitlog.txt %%fN%%_batch.log\n\n');
                if (run.floatBods(1).WamIGenMds == 0 && run.floatBods(1).Ngen > 0)
                    for n = 1:nbody
                        fprintf(fileID, 'del %%fN%%.pre %%fN%%.mod\n\n');
                    end
                    fprintf(fileID, 'set "t0=%%Time%%"\n');
                    fprintf(fileID, 'set "d0=%%Date%%"\n\n');
                    fprintf(fileID, '"%s\\wamit"\n\n', run.exePath);
                    fprintf(fileID, ':looppre\n');
                    fprintf(fileID, 'If not exist %%fN%%.pre then goto looppre\n\n');
                    fprintf(fileID, '"%s\\defmod"\n\n', run.exePath);
                    fprintf(fileID, ':loopmod\n');
                    fprintf(fileID, 'If not exist %%fN%%.mod then goto loopmod\n\n');
                    fprintf(fileID, '"%s\\wamit"\n\n', run.exePath);
                    fprintf(fileID, 'set "t1=%%Time%%"\n');
                    fprintf(fileID, 'set "d1=%%Date%%"\n\n');
                else
                    fprintf(fileID, 'set "t0=%%Time%%"\n');
                    fprintf(fileID, 'set "d0=%%Date%%"\n\n');
                    fprintf(fileID, '"%s\\wamit"\n\n', run.exePath);
                    fprintf(fileID, 'set "t1=%%Time%%"\n');
                    fprintf(fileID, 'set "d1=%%Date%%"\n\n');
                end
                
                fprintf(fileID, 'echo WAMIT Run: %%fN%% >> %%fN%%_batch.log\n');
                fprintf(fileID, 'echo Started: %%d0%% %%t0%% >> %%fN%%_batch.log\n');
                fprintf(fileID, 'echo Stopped: %%d1%% %%t1%% >> %%fN%%_batch.log\n');
                
                fclose(fileID);
            end
            
            % fnames
            filename = [run.folder '\fnames.wam'];
            fileID = fopen(filename, 'wt');
            
            fprintf(fileID, [run.runName '.cfg\n']);
            fprintf(fileID, [run.runName '.pot\n']);
            fprintf(fileID, [run.runName '.frc\n']);
            for n = 1:length(run.geoFiles)
                fprintf(fileID, [run.geoFiles{n} '.gdf\n']);
                if (run.CompDrift)
                    fprintf(fileID, [run.geoFiles{n} '.csf\n']);
                end
            end
            
            fclose(fileID);
            
            % Other files
            run.writeWamConfig;
            run.writeCfg;
            run.writePot;
            run.writeFrc;
            if (run.AutoCSF)
                for n = 1:nbody
                    if (isempty(run.floatBods(n).GeoFile))
                        geoFile = [run.runName num2str(n)];
                    else
                        geoFile = run.floatBods(n).GeoFile;
                    end
                    run.writeAutoCsf(run.floatBods(n), geoFile);
                end
            end
        end
        
        function [] = WriteRunIReadRAO(run, varargin)
            % WRITERUNIREADRAO This method is used if the WAMIT parameter,
            % IREADRAO is originally set to 1, so that the original run
            % just computes options 1-4, letting the user generate custom
            % RAOs before running again using this method, where IREADRAO
            % will be changed to a value of 2 or 3, in order to obtain
            % options 5-9.
            
            run.IReadRAO = 3; % Change so that WAMIT knows to read in the user-inputted RAOs:
                              % 2 - input real and imaginary parts - must
                              % be preceded by columns where the modulus
                              % and phase would be.
                              % 3 - provide just modulus and phase values.
                              % (see manual for more details)
            run.writeCfg; % Rewrite the .CFG file
            run.writeRAOs; % Write the user-generated RAO file for WAMIT to read.
            
        end % WriteRunIReadRAO
        
        function WriteRAOs(run)
            % WriteRAOs: Public version of writeRAOs - needed for creating
            % the drift force computation using the MoeSim class.
            run.writeRAOs;
        end
        
        function [] = Run(run, varargin)
            % Runs the batch file created to run Wamit with the system
            % command.
            % By default, Matlab waits for Wamit to finish running before
            % it continues.
            % 'Background' is an optional argument that runs Wamit in the
            % background so that Matlab is free to use.
            
            opts = checkOptions({'Background'}, varargin);
            
            if (opts)
                system(['cd "' run.folder '" && "' run.folder '\wam_run.bat" &']);
            else
                system(['cd "' run.folder '" && "' run.folder '\wam_run.bat"']);
            end
        end
    end
    
    methods (Access = protected)
        
        function [] = makeFolderAndSub(run, fold)
            % do nothing
        end
        
    end
    
    methods (Access = private)
        
        function [] = writeGdf(run, fb, geoFile, includeInt)
            
            geo = fb.PanelGeo;
            %%%%%% make a copy
            geo = PanelGeo(geo);
            geo.Translate(-fb.CenterRot); %translate to rotate about center of rotation.
            
            filename = [run.folder '\' geoFile '.gdf'];
            fileID = fopen(filename, 'wt');
            
            ulen = 1;
            g = 9.806650;
            
            Nwet = sum(geo.IsWets);
            Nint = sum(geo.IsInteriors);
            if (includeInt)
                count = Nwet;
            else
                count = Nwet - Nint;
            end
            
            fprintf(fileID, ['Model ' geoFile ', created: ' date '\n']);
            fprintf(fileID, '%8.4f %8.4f\n', ulen, g);
            fprintf(fileID, '%i %i \n', geo.Xsymmetry, geo.Ysymmetry);
            fprintf(fileID, '%i\n', count);
            
            pans = geo.Panels;
            
            for n = 1:geo.Count
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
                        fprintf(fileID, '\t%8.4f\t%8.4f\t%8.4f\n', verts(m,1), verts(m,2), verts(m,3));
                    end
                end
            end
            
            fclose(fileID);
            
            % write other files needed for the gdf (i.e. for new modes)
            if (~isempty(fb.WriteFileMeth))
                params = {fb.WriteParams, [fb.XYpos, fb.Zpos], fb.Angle};
                fb.WriteFileMeth(run.folder, geoFile, params);
            end
        end
        
        function [] = writeAutoCsf(run, fb, geoFile)
            % WRITEAUTOCSF Writes the .CSF file for automatic creation of
            % the control surface for computing drift forces.
            
            %%% NOTE: not sure that this will set the ocrrect panel size if
            %%% using the higher order method.. should check this later if
            %%% important (noted on 25/3/20)
            
            if isempty(run.boxCSF)
                run.boxCSF = false; % Use cylinder-shaped control surface by default
            end
            
            filename = [run.folder '\' geoFile '.csf']; % file must have same 
                                            % name as corresponding .GDF file.
            fileID = fopen(filename, 'wt');
            
            fprintf(fileID, ['Model ' geoFile ', created: ' date '\n']);
            fprintf(fileID, '1\n'); % ILOWHICSF (1 - higher order, 0 - lower order)
            fprintf(fileID, '%i %i \n',fb.PanelGeo.Xsymmetry,fb.PanelGeo.Ysymmetry); % Use same symmetry as .gdf file
            if isempty(fb.WamPanelSize) % i.e. if low order method was used for GDF
                csfPanelSize = fb.PanelGeo.AvgPanelArea*4; % use average panel size from GDF - this might be quite small so need to choose the multiplier here carefully
            else
                csfPanelSize = fb.WamPanelSize; % use same panel size as GDF
            end
            fprintf(fileID, '0 0 %4.2f\n',csfPanelSize); % NPATCSF (must be = 0 for auto),
            % ICDEF (must be = 0 for auto), PSZCSF (akin to PANEL_SIZE for the GDF)
            
            % Find required size of control surface (using PanelGeo if available)
            multip = 1.1; % multiplier used to ensure control surface encapsulates the body surface.
            
            if isempty(fb.PanelGeo) % create PanelGeo if does not exist, or derive verts if high order method used for gdf
                if run.FloatingBodies.WamILowHi == 0 % if low order method used for gdf
                    fb.PanelGeo = Wamit_readGdf(run.Folder, geoFile);
                    for i = 1:size(fb.PanelGeo.Panels,2)
                        panVerts_x(:,i) = fb.PanelGeo.Panels(i).Vertices(:,1);
                        panVerts_y(:,i) = fb.PanelGeo.Panels(i).Vertices(:,2);
                        panVerts_z(:,i) = fb.PanelGeo.Panels(i).Vertices(:,3);
                    end
                elseif run.FloatingBodies.WamILowHi == 1 % if high order used for gdf
                    verts = Wamit_readGdfHi(run.Folder, geoFile);
                    for i = 1:size(verts,1)
                        vertsSize(i) = size(verts{i,1},1);
                    end
                    panVerts_x = zeros(max(vertsSize),size(verts,1));
                    for i = 1:size(verts,1)
                        panVerts_x(1:length(verts{i,1}),i) = verts{i,1}(:,1);
                        panVerts_y(1:length(verts{i,1}),i) = verts{i,1}(:,2);
                        panVerts_z(1:length(verts{i,1}),i) = verts{i,1}(:,3);
                    end
                end
            else % i.e. if the panelGeo is already defined
                for i = 1:size(fb.PanelGeo.Panels,2)
                    panVerts_x(:,i) = fb.PanelGeo.Panels(i).Vertices(:,1);
                    panVerts_y(:,i) = fb.PanelGeo.Panels(i).Vertices(:,2);
                    panVerts_z(:,i) = fb.PanelGeo.Panels(i).Vertices(:,3);
                end
            end
            
            x_max = max(abs(max(max(panVerts_x))),abs(min(min(panVerts_x))));
            y_max = max(abs(max(max(panVerts_y))),abs(min(min(panVerts_y))));
            rad = max(x_max,y_max)*multip;
            
            z_max = min(min(panVerts_z)); % I think this should be the minimum z-coordinate
            % as measured from the free
            % surface.
            dep = abs(z_max)*multip;
            npart = 0; % must = 0 for cylinder by WAMIT convention (see section 11.5 of v7.3 manual)
            
            if run.boxCSF 
                rad = 0; % For box, need to set rad = 0 since vertices below define the x-y proections of surface
                npart = 1; % x-y vertices definition to follow below...
            end
            
            fprintf(fileID, '%4.2f %4.2f\n',rad,dep); % RADIUS (>0 for cylinder, <=0 for box, 
                                        % box needs additional info, see wamit manual), 
                                        % DEPTH (vertical coverage of control surface)
            fprintf(fileID, '%i\n', npart); % NPART (=0 if using cylinder and only exists one waterline, =1 for a single box-shaped surface.)
            
            if run.boxCSF
                % Find x-y vertices of control surface
                xu = max(max(panVerts_x))*multip; xl = min(min(panVerts_x))*multip;
                yu = max(max(panVerts_y))*multip; yl = min(min(panVerts_y))*multip;
                
                fprintf(fileID, '4\n'); % Currently set up to allow only 4 vertices to form box around body surface
                
                fprintf(fileID, '%4.4f  %4.4f\n', xu, yu); % Coordinates must proceed in a counter-clockwise
                fprintf(fileID, '%4.4f  %4.4f\n', xl, yu); % direction around the edge of the control surface,
                fprintf(fileID, '%4.4f  %4.4f\n', xl, yl); % as viewed from above.
                fprintf(fileID, '%4.4f  %4.4f\n', xu, yl);
    
            end
            
            fclose(fileID);
            
        end % writeAutoCsf
        
        function [] = writeWamConfig(run)
            filename = [run.folder '\config.wam'];
            fileID = fopen(filename, 'wt');

            fprintf(fileID, '! generic configuration file:  config.wam\n');
            fprintf(fileID, 'RAMGBMAX=%4.1f\n', run.ramGB);
            fprintf(fileID, 'NCPU=%i\n', run.ncpu);
            fprintf(fileID, 'USERID_PATH = %s\n', run.useridPath);
            fprintf(fileID, 'LICENSE_PATH=\\wamitv7\\license\n');
            fclose(fileID);
        end
        
        function [] = writeCfg(run)
            Nbod = length(run.floatBods);
            filename = [run.folder '\' run.runName '.cfg'];
            fileID = fopen(filename, 'wt');

            fprintf(fileID, 'IDELFILES = 2\n');
            fprintf(fileID, 'IALTPOT = 2\n');
            if (run.compDrift & ismember(1,run.driftOption))
            fprintf(fileID, 'IALTCSF = 2\n'); % Alternative 2 must be used 
                         % if controls surfaces are very close to the body 
                         % surface. (it includes the waterline integral)
            end
            if (run.floatBods(1).ISurfPan)
                irr = 1;
            else
                irr = 0;
            end
            if (run.floatBods(1).MakeSurfPan)
                irr = 3;
            end
            fprintf(fileID, 'IRR = %i\n', irr);
            fprintf(fileID, 'ILOG = 1\n');
            fprintf(fileID, 'IALTFRC = 2\n');
            if (~isempty(run.fieldArray))
                fprintf(fileID, 'IFIELD_ARRAYS = 1\n');
            end
            fprintf(fileID, 'IFORCE = 1\n');
            fprintf(fileID, 'ILOWHI = %i\n', run.floatBods(1).WamILowHi);
            if (run.floatBods(1).WamILowHi)
                if isempty(run.floatBods(1).WamPanelSize)
                    error('WamPanelSize not set for WamILowHi = 1');
                end
                fprintf(fileID, 'PANEL_SIZE = %i\n', run.floatBods(1).WamPanelSize); 
            end
            if (run.floatBods(1).WamILowHi)
                fprintf(fileID, 'ILOWGDF = 0\n');
                if (sum((run.floatBods(1).WamSpline)~=0))
                    fprintf(fileID, 'KSPLIN = %i\n', run.floatBods(1).WamSpline(1));
                    fprintf(fileID, 'IQUADO = %i\n', run.floatBods(1).WamSpline(2));
                    fprintf(fileID, 'IQUADI = %i\n', run.floatBods(1).WamSpline(3));
                end
                    
            end
            if run.floatBods(1).SurfAboveZ0
                fprintf(fileID, 'ITRIMWL = %i\n', 10);
%                 fprintf(fileID, 'ITRIMWL = %i\n', Nbod);
                if Nbod > 1
                    for n = 1:Nbod
                        fprintf(fileID, 'XTRIM(%i) = 0.0 0.0 0.0\n', n);
                    end
                else
                    fprintf(fileID, 'XTRIM = 0.0 0.0 0.0\n');
                end
            end
            if (~all([isempty(run.fieldArray) isempty(run.fieldPoints) isempty(run.cylArray)]))
                fprintf(fileID, 'INUMOPT6 = 1 \n');
            end
            if run.computeBody && isempty(run.fieldPointAsBodyDist)
                fprintf(fileID, 'INUMOPT5 = 1 \n');
                if (run.floatBods(1).WamILowHi > 0) || ~isempty(run.floatBods(1).CompGeoFile)
                    ipnlbpt = -4;
                else
                    ipnlbpt = 0;
                end
                fprintf(fileID, 'IPNLBPT = %i \n', ipnlbpt);
            end
            % Generalized modes
            if Nbod > 1
                for n = 1:Nbod
                    fprintf(fileID, 'NEWMDS(%i) = %i\n', n, run.floatBods(n).Ngen);
                    fprintf(fileID, 'IGENMDS(%i) = %i\n', n, run.floatBods(n).WamIGenMds);
                end
            else
                fprintf(fileID, 'NEWMDS = %i\n', run.floatBods(1).Ngen);
                fprintf(fileID, 'IGENMDS = %i\n', run.floatBods(1).WamIGenMds);
            end
            
            
            
            dipoles = run.floatBods.WamDipoles;
            if ~isempty(dipoles)
                if iscell(dipoles)
                    % New method requires cell array of panel indices
                    for n = 1:Nbod
                        dipoles = run.floatBods(n).WamDipoles;
                        if ~isempty(dipoles)
                            if Nbod == 1
                                fprintf(fileID, 'NPDIPOLE =');
                            else
                                fprintf(fileID, 'NPDIPOLE(%i) =', n);
                            end
                            for m = 1:length(dipoles)
                                dipm = dipoles{m};
                                if length(dipm) == 0
                                    % skip
                                elseif length(dipm) == 2
                                    fprintf(fileID, ' (%i %i)', dipm(1), dipm(2));
                                elseif length(dipm) == 1
                                    fprintf(fileID, ' %i', dipm(1));
                                else
                                    error('Unrecognized array size in Wamit Dipole panels');
                                end
                            end
                            fprintf(fileID, '\n');
                        end
                    end
                else
                    % Old method: Doesn't use cell array for capturing
                    % dipoles. Can't do two individual panels - i.e.
                    % assumes [10 12] is the range 10-12, rather than
                    % individual panels 10 and 12. Also, can't do
                    % combinations of ranges, and individual panels.
                    for n = 1:Nbod
                        dipoles = run.floatBods(n).WamDipoles;
                        if ~isempty(dipoles)
                            if length(dipoles) == 2
                                if Nbod == 1
                                    fprintf(fileID, 'NPDIPOLE = (%i %i)\n', dipoles(1), dipoles(2));
                                else
                                    fprintf(fileID, 'NPDIPOLE(%i) = (%i %i)\n', n, dipoles(1), dipoles(2));
                                end
                            else
                                if Nbod == 1
                                    fprintf(fileID, 'NPDIPOLE =');
                                else
                                    fprintf(fileID, 'NPDIPOLE(%i) =', n);
                                end
                                for m = 1:length(dipoles)
                                    fprintf(fileID, ' %i', dipoles(m));
                                end
                                fprintf(fileID, '\n');
                            end
                        end
                    end
                end
            end
            if run.iReadRAO > 1 % i.e. if running WAMIT for the second time to evaluate options 5-9
                ipoten = 0; % Don't solve for potentials
            else
                ipoten = 1; % Solve for potentials
            end
            fprintf(fileID, 'IPOTEN = %i\n', ipoten);
            fprintf(fileID, 'IREADRAO = %i\n', run.iReadRAO); % This will = 0 if user has not switched it to a value >0.
            fprintf(fileID, 'ISCATT = 0\n');
            if (run.useDirect)
                % User input to use the direct solver
                fprintf(fileID, 'ISOLVE = 1\n');
            elseif ~isempty(run.blockItCnt)
                fprintf(fileID, 'ISOLVE = %i\n', run.blockItCnt);
            else
                if (run.floatBods(1).WamILowHi)
                    % Use direct solver for higher order panel method
                    fprintf(fileID, 'ISOLVE = 1\n');
                else
                    % Use iterative solver for lower order panel method
                    fprintf(fileID, 'ISOLVE = 0\n');
                end
            end
            if (run.computeBody && run.computeVelocity && ~run.floatBods(1).WamILowHi) || ...
                    (~run.floatBods(1).WamILowHi && run.CompDrift)
                fprintf(fileID, 'ISOR = 1\n');
            else
                fprintf(fileID, 'ISOR = 0\n');
            end
            fprintf(fileID, 'MAXITT = %i\n', run.maxItt);
            fprintf(fileID, 'MAXMIT = 8\n');
            fprintf(fileID, 'MONITR = 0\n');
            fprintf(fileID, 'NOOUT = 1  1  1  1  0  0  0  0  0\n');
            fprintf(fileID, 'NUMHDR = 1\n');
            fprintf(fileID, 'IPLTDAT = 1\n');
            fprintf(fileID, 'NUMNAM = 0\n');
            fprintf(fileID, 'SCRATCH_PATH = %s\n', run.scratchPath);

            fclose(fileID);
        end
        
        function [] = writeFrc(run)
            filename = [run.folder '\' run.runName '.frc'];
            fileID = fopen(filename, 'wt');
            
            % header
            fprintf(fileID, ['.frc file   Series: '  '    Run: ' run.runName '\n']);
            % options
            ibp = 1;
            ifldp = 1;
            if (run.computeVelocity)
                ifldp = 3;
                ibp = 3;
            end
            if ~run.computeBody
                ibp = 0;
            elseif ~isempty(run.fieldPointAsBodyDist)
                ibp = 0;
            end
            if (all([isempty(run.fieldArray) isempty(run.fieldPoints) isempty(run.cylArray)]))
                ifldp = 0;
            end
            % 1 - added mass and damping
            % 2 - haskind exciting forces
            % 3 - diffraction potential exciting forces
            % 4 - body motions
            % 5 - pressure and velocity on body
            % 6 - pressure velocity: 0 - no output; 1/-1 - pressure
            % potential/source; 2/-2 - velocity potential/source; 3/-3 -
            % both potential/source
            % 7 - mean drift force using control surface
            % 8 - mean drift forces using conservation of momentum 
            % 9 - mean drift forces using pressure integration
            iradf = 1;
            if (run.computeFK)
                fk = 2;
            else
                fk = 1;
            end
            if (run.useHaskind)
                iexH = fk;
                iexD = 0;
            else
                iexH = 0;
                iexD = fk;
            end
            if (~run.solveRad)
                iradf = 0;
                iexH = 0;
                iexD = 1;
            end
            imot = 2; % Output body motions (& use these for computation of drift forces (and vels and pressures if applicable)) Haskind - 1, diffraction - 2
            %%% MAY WANT TO CHANGE THIS LATER TO TURN THIS OPTION OFF WHEN
            %%% IREADRAO = 2 or 3 FOR DRIFT FORCE CALCULATIONS.
            icompDrift = [0 0 0]; % Initialise vector of drift options
            if (run.compDrift) % Using control surface
                for i = 1:length(run.driftOption)
                    icompDrift(run.driftOption(i)) = 1; % (switch each driftOption on)
                end
            end
            fprintf(fileID, '%i  %i  %i  %i  %i  %i  %i  %i  %i\n', iradf, iexH, iexD, imot, ibp, ifldp, icompDrift);  
            % rho
            fprintf(fileID, '%8.4f\n', run.rho);
            % cg
            nbody = length(run.floatBods);
            ndof = 0;
            for n = 1:nbody
                fb = run.floatBods(n);
                cg = fb.Cg;
                cR = fb.CenterRot;
                cg = cg - cR;
                fprintf(fileID, '%8.4f\t%8.4f\t%8.4f\t', cg(1), cg(2), cg(3));
                mod = fb.Modes;
                ndof = ndof + 6 + fb.Ngen;
            end
            % Get Mass, Damping and Stiffness matrices from floating bodies.
            % Can't do coupling between floating bodies yet...
            M = zeros(ndof, ndof);
            D = zeros(ndof, ndof);
            K = zeros(ndof, ndof);
            istart = 1;
            for n = 1:nbody
                fb = run.floatBods(n);
                istop = istart + 5 + fb.Ngen;
                M(istart:istop, istart:istop) = fb.M;
                D(istart:istop, istart:istop) = fb.Dpto;
                K(istart:istop, istart:istop) = fb.K;
                istart = istop + 1;
            end
            % IMASS - include mass matrix (1) or not (0)?  
            fprintf(fileID, '\n1\n');
            % Mass matrix
            for n = 1:ndof
                for m = 1:ndof
                    fprintf(fileID, '%8.4f\t', M(n,m)); 
                end
                fprintf(fileID, '\n');
            end
            % IDAMP - include damping matrix (1) or not (0)?  
            fprintf(fileID, '1\n');
            % Damping matrix
            for n = 1:ndof
                for m = 1:ndof
                    fprintf(fileID, '%8.4f\t', D(n,m)); 
                end
                fprintf(fileID, '\n');
            end
            % ISTIF - include stiffness matrix (1) or not (0)?  
            fprintf(fileID, '1\n');
            % stiffness matrix
            for n = 1:ndof
                for m = 1:ndof
                    fprintf(fileID, '%8.4f\t', K(n,m)); 
                end
                fprintf(fileID, '\n');
            end
            % NBETAH - number of Haskind wave headings. Haskind wave headings are in
            % addition to diffraction force headings in the potential control file.
            % Haskind diffraction forces only require radiation potential.
            fprintf(fileID, '0\n');
      
            %{
            % Field points Method for 6.0
            points = run.fieldPoints;
            farray = run.fieldArray;
            if (~isempty(farray))
                pointsF = points;
                NF = size(pointsF, 1);
                
                pointsA = farray.GetPointsList;
                NA = size(pointsA, 1);
                
                points = zeros(NF + NA, 3);
                if (NF > 0)
                    points(1:NF,:) = pointsF;
                end
                points(NF+1:end,:) = pointsA;
            end
            N = size(points, 1);
            fprintf(fileID, '%i\n', N);
            % explicit field points
            for n = 1:N
                fprintf(fileID, '%8.4f\t%8.4f\t%8.4f\t\n', points(n,1), points(n,2), points(n,3));
            end
            %}
            
            % Field points Method for 6.4 and above...
            % NFIELD - number of explicitly specified field points
            if (~isempty(run.cylArray))
                if (~isempty(run.fieldPoints))
                    error('Cannot create WAMIT run with cylinder array and field points');
                end
                points = run.cylArray.GetPoints(run.h);
            else
                points = run.fieldPoints;
            end
            N = size(points, 1);
            fprintf(fileID, '%i\n', N);
            % explicit field points
            for n = 1:N
                fprintf(fileID, '%8.4f\t%8.4f\t%8.4f\t\n', points(n,1), points(n,2), points(n,3));
            end
            % NFIELD_ARRAYS - number of seperate arrays. In cfg file, IFIELD_ARRAY = 1.
            farray = run.fieldArray;
            if (~isempty(farray))
                fprintf(fileID, '%i\n', 1);
                % ITANKFLD - outside tank (0), inside tank (1)
                fprintf(fileID, '0\n'); 
                np = farray.NumberPoints;
                start = farray.Start;
                deltas = farray.Deltas;
                % x-direction: number of points, start point, delta
                fprintf(fileID, '%i\t%8.4f\t%8.4f\n', np(1), start(1), deltas(1));
                % y-direction: number of points, start point, delta
                fprintf(fileID, '%i\t%8.4f\t%8.4f\n', np(2), start(2), deltas(2));
                % z-direction: number of points, start point, delta
                fprintf(fileID, '%i\t%8.4f\t%8.4f\n', np(3), start(3), deltas(3));
            end
            
            fclose(fileID);
        end
        
        function [] = writeRAOs(run)
            % WRITERAOS Writes RAOs computed outside of WAMIT to a .RAO
            % file of the same name as the FRC file. Used for input to
            % WAMIT in order to compute options 5-9 (pressures, velocities 
            % or drift forces).
            
            % First, read the mode indices so know which DoFs are which in
            % the motionRAOs matrix
            try
                file_data = importdata([run.folder '/' run.runName '.4']);
            catch
                error('WamitRunCondition:writeRAOs:noRAOdemoFile',...
                    'The demo RAO file does not exist.');
            end
            mode_indices = unique(file_data.data(:,3));
            if length(mode_indices)~=size(run.motionRAOs,3)
                error('WamitRunCondition:writeRAOs:motionIndicesDoNotMatch',['The mode indices from the original .4 file do ' ... 
                    'not match those given in the motionRAOs structure.']);
            end
%             mode_indices = mode_indices(mode_indices<=12); % DoFs 1-12 as numbered by wamit will always refer to the two hulls,
                                                           % any additional indices will refer to damping plates.  
            
            % Now, open and write the RAO file to be read in by WAMIT            
            filename = [run.folder '\' run.runName '.rao'];
            fileID = fopen(filename, 'wt');
            
            % header
            fprintf(fileID, ['.rao file   Series: '  '    Run: ' run.runName '\n']);
            
            % print RAOs
            for i = 1:length(run.t)
                for j = 1:length(run.beta)
                    for k = 1:length(mode_indices)
                        fprintf(fileID, '%8.4f\t%8.4f\t%i\t%8.6e\t%8.6e\n', ...
                            run.t(i),run.beta(j)*180/pi,mode_indices(k),...
                            abs(run.motionRAOs(i,j,k)),angle(run.motionRAOs(i,j,k))*180/pi);
                        % period, angle (in degrees), DoF, abs(RAO), angle(RAO) in degrees
                    end
                end
            end
            fclose(fileID);
            
        end % writeRAOs
        
        function [] = writePot(run)
            filename = [run.folder '\' run.runName '.pot'];
            fileID = fopen(filename, 'wt');
            
            % header
            fprintf(fileID, ['.pot file   Series: '  '    Run: ' run.runName '\n']);
            % HBOT - dimensional water depth (infinite depth: -1)
            h_ = -1;
            if (~isinf(run.h))
                h_ = run.h;
            end
            fprintf(fileID, '%8.4f\n', h_);
            % IRAD IDIFF (1 = solve all 6 dof, 0 = solve only specified modes, -1 =
            % don't solve)
            irad = -1;
            if (run.solveRad)
                irad = 0;
            end
            idiff = -1;
            if (run.solveDiff)
                idiff = 0;
            end
            fprintf(fileID, '%i           %i\n', irad, idiff);
            % NPER - number of periods
            t_ = run.t;
            nT = length(t_);
            nTw = nT;
            if (run.incInfFreq)
                nTw = nTw + 1;
            end
            if (run.incZeroFreq)
                nTw = nTw + 1;
            end
            fprintf(fileID, '%i		Number of periods to be analyzed\n', nTw);
            if (run.incZeroFreq)
                fprintf(fileID, '%8.4f\n', -1);
            end
            if (run.incInfFreq)
                fprintf(fileID, '%8.4f\n', 0);
            end
            % periods
            for n = 1:nT
                fprintf(fileID, '%8.4f\n', t_(n));
            end
            % NBETA - number of directions
            bet = run.beta;
            nB = length(bet);
            fprintf(fileID, '%i		Number of direction headdings to be analyzed\n', nB);
            % directions
            for n = 1:nB
                fprintf(fileID, '%8.4f\n', round(1e4*180/pi*bet(n))/1e4);
            end
            % NBODY - number of bodies
            nbody = length(run.floatBods);
            fprintf(fileID, '%i\n', nbody);
            for n = 1:nbody
                % Name of gdf
                fb = run.floatBods(n);
                fprintf(fileID, [run.geoFiles{n} '.gdf\n']);
                % x, y, z, and angle (degrees of body x-axis relative to global x-axis) of
                % body n
                %%%%
                cR = fb.CenterRot;
                pos = [fb.XYpos(1) fb.XYpos(2) fb.Zpos];
                pos = pos + cR; % add cr bc body moved so cr is the origin. 
                fprintf(fileID, '%8.4f\t%8.4f\t%8.4f\t%8.4f\n', pos(1), pos(2), pos(3), fb.Angle);
                % modes to be computed of body n
                modes = fb.Modes;
                vect = modes.Vector;
                fprintf(fileID, '%i %i %i %i %i %i\n', vect(1:6));
            end

            fclose(fileID);
        end
    end
    
    methods (Static)
        function [] = MakeSuperBatch(runs)
            N = length(runs);
            
            for n = 1:N
                if (~isa(runs(n), 'Wamit_RunCondition'))
                    error('All runs must be of type Wamit_RunCondition');
                end
            end
            
            path = runs(1).Folder;
            
            for n = 2:N
                if (~strcmp(path, runs(n).Folder))
                    error('All runs must have the same Folder');
                end
            end
            
            filename = [path '\wam_superRun.bat'];
            fileID = fopen(filename, 'wt');

            fprintf(fileID, 'set "t0=%%Time%%"\n');
            fprintf(fileID, 'set "d0=%%Date%%"\n\n');
            
            for n = 1:N
                fprintf(fileID, ['cd '  runs(n).RunName '\n']);
                fprintf(fileID, 'wam_run.bat\n\n');        
                fprintf(fileID, [':loop' num2str(n) '\n']);
                fprintf(fileID, ['If not exist ' runs(n).RunName '_batch.log then goto loop' num2str(n) '\n\n']);
                fprintf(fileID, 'cd ..\n\n');
            end
            
            fprintf(fileID, 'set "t1=%%Time%%"\n');
            fprintf(fileID, 'set "d1=%%Date%%"\n\n');
            
            fprintf(fileID, 'echo WAMIT Super Run: %%fN%% >> ws_batch.log\n');
            fprintf(fileID, 'echo Started: %%d0%% %%t0%% >> ws_batch.log\n');
            fprintf(fileID, 'echo Stopped: %%d1%% %%t1%% >> ws_batch.log\n');

            fclose(fileID);
        end
    end
    
end