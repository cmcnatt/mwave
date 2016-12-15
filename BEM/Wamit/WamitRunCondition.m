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
        scratchPath;
        useridPath;
        ncpu;
        ramGB;
        maxItt;
        useDirect;
        useHaskind;
        computeFK;
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
        ComputeVelocity;    % Indicates whether velocity should be computed at field points
        WamitPath;          % Path location of wamit.exe
        ScratchPath;        % Path location of wamit scratch folder
        UseridPath;         % Path location of userid.wam license file
        NCPU;               % The number of cpus to be used in the Wamit computation
        RAMGBmax;           % The max RAM to be used in the Wamit computation
        MaxItt;
        UseDirectSolver;
        UseHaskind;
        ComputeFK;
    end

    methods
        
        function [run] = WamitRunCondition(folder, runName)
            run.rho = 1000;
            run.computeVelocity = false;
            run.computeBody = false;
            run.solveDiff = true;
            run.solveRad = true;
            run.exePath = 'C:\wamitv7';
            run.scratchPath = '\wamitv7\scratch';
            run.useridPath = '\wamitv7';
            run.ncpu = 1;
            run.ramGB = 2;
            run.maxItt = 35;
            run.useDirect = false;
            run.useHaskind = false;
            run.computeFK = false;
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
            if ((length(val) > 1) && run.computeBody)
                error('Cannot compute body points for an array of floating bodies');
            end
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
        
        function [] = WriteRun(run, varargin)
            % Writes the Input files (except for the .gdf) for a WAMIT run
            
            [opts] = checkOptions({'NoBatch'}, varargin);
            noBat = opts(1);
            
            if (exist(run.folder, 'file') ~= 7)
                error('Designated run folder does not exist');
            end
                        
            % .gdf
            nbody = length(run.floatBods);
            run.geoFiles = cell(1, nbody);
            for n = 1:nbody
                if (isempty(run.floatBods(n).GeoFile))
                    geoFile = [run.runName num2str(n)];
                    run.writeGdf(run.floatBods(n), geoFile, run.floatBods(n).ISurfPan);
                    run.geoFiles{n} = geoFile;
                    %run.floatBods(n).WriteGdf(run.folder, fileName);
                else
                    run.geoFiles{n} = run.floatBods(n).GeoFile;
                end
            end
            
            if (~noBat)
                filename = [run.folder '\wam_run.bat'];
                fileID = fopen(filename, 'wt');
                
                fprintf(fileID, ['set "fN=' run.runName '"\n\n']);
                fprintf(fileID, 'del %%fN%%.out %%fN%%.1 %%fN%%.2 %%fN%%.3 %%fN%%.4 ');
                fprintf(fileID, '%%fN%%.5 %%fN%%.6p %%fN%%.6vx %%fN%%.6vy %%fN%%.6vz ');
                fprintf(fileID, '%%fN%%.fpt %%fN%%.p2f errorf.log ');
                fprintf(fileID, 'errorp.log rgkinit.txt rgklog.txt wamitlog.txt %%fN%%_batch.log\n\n');
                if (run.floatBods(1).WamIGenMds == 0 && run.floatBods(1).Ngen > 0)
                    for n = 1:nbody
                        fprintf(fileID, 'del %%fN%%.pre %%fN%%.mod\n\n');
                    end
                    fprintf(fileID, 'set "t0=%%Time%%"\n');
                    fprintf(fileID, 'set "d0=%%Date%%"\n\n');
                    fprintf(fileID, '%s\\wamit\n\n', run.exePath);
                    fprintf(fileID, ':looppre\n');
                    fprintf(fileID, 'If not exist %%fN%%.pre then goto looppre\n\n');
                    fprintf(fileID, '%s\\defmod\n\n', run.exePath);
                    fprintf(fileID, ':loopmod\n');
                    fprintf(fileID, 'If not exist %%fN%%.mod then goto loopmod\n\n');
                    fprintf(fileID, '%s\\wamit\n\n', run.exePath);
                    fprintf(fileID, 'set "t1=%%Time%%"\n');
                    fprintf(fileID, 'set "d1=%%Date%%"\n\n');
                else
                    fprintf(fileID, 'set "t0=%%Time%%"\n');
                    fprintf(fileID, 'set "d0=%%Date%%"\n\n');
                    fprintf(fileID, '%s\\wamit\n\n', run.exePath);
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
            end
            
            fclose(fileID);
                        
            run.writeWamConfig;
            run.writeCfg;
            run.writePot;
            run.writeFrc;
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
                system(['cd ' run.folder ' && ' run.folder '\wam_run.bat &']);
            else
                system(['cd ' run.folder ' && ' run.folder '\wam_run.bat']);
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
            filename = [run.folder '\' run.runName '.cfg'];
            fileID = fopen(filename, 'wt');

            fprintf(fileID, 'IDELFILES = 2\n');
            fprintf(fileID, 'IALTPOT = 2\n');
            if (run.floatBods(1).ISurfPan)
                irr = 1;
            else
                irr = 0;
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
                fprintf(fileID, 'PANEL_SIZE = %i\n', run.floatBods(1).WamPanelSize); 
            end
            if (run.floatBods(1).WamILowHi)
                fprintf(fileID, 'ILOWGDF = 0\n');
            end
            if (~all([isempty(run.fieldArray) isempty(run.fieldPoints) isempty(run.cylArray)]))
                fprintf(fileID, 'INUMOPT6 = 1 \n');
            end
            if (run.computeBody)
                fprintf(fileID, 'INUMOPT5 = 1 \n');
                if run.floatBods(1).WamILowHi > 0
                    ipnlbpt = -4;
                else
                    ipnlbpt = 0;
                end
                fprintf(fileID, 'IPNLBPT = %i \n', ipnlbpt);
            end
            % Generalized modes
            Nbod = length(run.floatBods);
            if (Nbod > 1)
                for n = 1:Nbod
                    fprintf(fileID, 'NEWMDS(%i) = %i\n', n, run.floatBods(n).Ngen);
                    fprintf(fileID, 'IGENMDS(%i) = %i\n', n, run.floatBods(n).WamIGenMds);
                end
            else
                fprintf(fileID, 'NEWMDS = %i\n', run.floatBods(1).Ngen);
                fprintf(fileID, 'IGENMDS = %i\n', run.floatBods(1).WamIGenMds);
            end
            % Dipole patches
            % NPDIPOLE(2)=(5-8)
            if Nbod > 1
                for n = 1:Nbod
                    dipoles = run.floatBods(n).WamDipoles;
                    if ~isempty(dipoles)
                        fprintf(fileID, 'NPDIPOLE(%i) = ', n);
                        for m = 1:length(dipoles)
                            fprintf(fileID, '%i', dipoles(m));
                            if m ~= length(dipoles)
                                fprintf(fileID, ', ');
                            end
                        end
                        fprintf(fileID, '\n');
                    end
                end
            else
                dipoles = run.floatBods.WamDipoles;
                if ~isempty(dipoles)
                    fprintf(fileID, 'NPDIPOLE = ');
                    for m = 1:length(dipoles)
                        fprintf(fileID, '%i', dipoles(m));
                        if m ~= length(dipoles)
                            fprintf(fileID, ', ');
                        end
                    end
                    fprintf(fileID, '\n');
                end
            end
            fprintf(fileID, 'IPOTEN = 1\n');
            fprintf(fileID, 'ISCATT = 0\n');
            if (run.useDirect)
                % User input to use the direct solver
                fprintf(fileID, 'ISOLVE = 1\n');
            else
                if (run.floatBods(1).WamILowHi)
                    % Use direct solver for higher order panel method
                    fprintf(fileID, 'ISOLVE = 1\n');
                else
                    % Use iterative solver for lower order panel method
                    fprintf(fileID, 'ISOLVE = 0\n');
                end
            end
            if (run.computeBody && run.computeVelocity)
                fprintf(fileID, 'ISOR = 1\n');
            else
                fprintf(fileID, 'ISOR = 0\n');
            end
            fprintf(fileID, 'MAXITT = %i\n', run.maxItt);
            fprintf(fileID, 'MAXMIT = 8\n');
            fprintf(fileID, 'MONITR = 0\n');
            fprintf(fileID, 'NOOUT = 1  1  1  1  0  0  0  0  0\n');
            fprintf(fileID, 'NUMHDR = 1\n');
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
            if (~run.computeBody)
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
            % both potenital/source
            % 7 - mean drift forces
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
            fprintf(fileID, '%i  %i  %i  0  %i  %i  0 0  0\n', iradf, iexH, iexD, ibp, ifldp);  
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
                idiff = 1;
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
                %fprintf(fileID, '%8.4f\t%8.4f\t%8.4f\t%8.4f\n', fb.XYpos(1), fb.XYpos(2), fb.Zpos, fb.Angle);
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