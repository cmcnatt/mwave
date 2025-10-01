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
classdef NemohRunCondition < IBemRunCondition
    % Creates input files for a Nemoh run
    
    properties (Access = private)
        nx;
        nNode;
        nf;
        nPan;
        omegaLim;
        nOmega;
        betaLim;
        nBeta;
        computeHS;
        useCylSurf;        
    end

    properties (Dependent)
        OmegaLims;          % Lower and upper limit of the radial frequencies
        OmegaCount;         % Number of radial frequencies
        BetaLims;           % Lower and upper limit of the wave direction (radians)
        BetaCount;          % Number of wave directions
        T;                  % The wave periods (can't set)
        Beta;               % The wave directions (can't set)
        FloatingBodies;     % FloatingBodies (plural because of multiple bodies)
        ComputeHS;          % Indicates whether to use the Nemoh Mesh.exe to compute the hydrostatics
    end
    
    methods
        
        function [run] = NemohRunCondition(folder)
            run.exePath = ['/Users/rcoe/Nemoh/bin'];
            run.rho = 1000;
            run.nOmega = 0;
            run.nBeta = 0;
            run.geoFiles = [];
            run.computeHS = true;
            run.fieldArray = [];
            run.cylArray = [];
            run.useCylSurf = true;
            if (nargin == 0)
                run.folder = '.';
            else 
                run.folder = folder;
                run.makeFolderAndSub(folder);
            end
        end
                
        function [ol] = get.OmegaLims(run)
            % Get the lower and upper freguency limits
            ol = run.omegaLim;
        end
        function [run] = set.OmegaLims(run, ol)
            % Get the lower and upper freguency limits
            
            if (length(ol) ~= 2)
                error('The limits must be a vector of length 2');
            end
            
            if (ol(2) < ol(1))
                error('The first element is the lower limit, which must be less than the upper limit');
            end
            
            if any(ol <= 0)
                error('The limits must be positive');
            end
            
            run.omegaLim = ol;
        end
        
        function [oc] = get.OmegaCount(run)
            % Get the number of frequencies
            oc = run.nOmega;
        end
        function [run] = set.OmegaCount(run, oc)
            % Set the number of frequencies
            if (oc <= 0)
                error('The number of frequnecies must be a positive integer');
            end
            
            if (~isInt(oc))
                error('The number of frequnecies must be a positive integer');
            end
            
            run.nOmega = oc;
        end
        
        function [t] = get.T(run)
            % Get the wave periods
            t = 2*pi./linspace(run.omegaLim(1), run.omegaLim(2), run.nOmega);
        end
        
        function [bl] = get.BetaLims(run)
            % Set the wave headings to be run (radians)           
            bl = run.betaLim;
        end
        function [run] = set.BetaLims(run, bl)
            % Set the wave headings to be run (radians)
            
            if (length(bl) ~= 2)
                error('The limits must be a vector of length 2');
            end
            
            if (bl(2) < bl(1))
                error('The first element is the lower limit, which must be less than the upper limit');
            end
            
            if any(bl < 0)
                error('The limits must be positive');
            end
            
            run.betaLim = bl;
        end
        
        function [bc] = get.BetaCount(run)
            % Get the number of incident directions
            bc = run.nBeta;
        end
        function [run] = set.BetaCount(run, bc)
            % Set the number of incident directions
            if (bc <= 0)
                error('The number of directions must be a positive integer');
            end
            
            if (~isInt(bc))
                error('The number of directions must be a positive integer');
            end
            
            run.nBeta = bc;
        end
        
        function [beta] = get.Beta(run)
            % Get the wave directions
            beta = linspace(run.betaLim(1), run.betaLim(2), run.nBeta);
        end
        
        function [fbs] = get.FloatingBodies(run)
            % Get the array of floating bodies
            fbs = run.floatBods;
        end
        function [run] = set.FloatingBodies(run, fbs)
            % Set the array of floating bodies
            pass = 1;
            
            for n = 1:length(fbs)
                if (~isa(fbs(n), 'FloatingBody'))
                    pass = 0;
                end
                if (fbs(n).Ngen > 0)
                    error('Generalized modes not yet implemented with Nemoh');
                end
            end
            if (pass)
                run.floatBods = fbs;
            else
                error('All Floating Bodies must be of type FloatingBody');
            end    
        end
        
        function [hs] = get.ComputeHS(run)
            % Indicates whether or not to compute HS with Mesh.exe
            hs = run.computeHS;
        end
        function [run] = set.ComputeHS(run, hs)
            % Indicates whether or not to compute HS with Mesh.exe
            if (~isBool(hs))
                error('The hydrostatic indicator must be a boolean');
            end
            
            run.computeHS = hs;
        end
                        
        function [] = WriteRun(run)
            % Writes the files for the Nemoh run
                        
                        
            % write Mesh cal - just use this to compute Hydrostatic matrix
            % for now, only works on first body
            if (run.computeHS)
                run.writeMeshCalAndFile(run.floatBods(1));
            end
            
            % write the Nemoh Cal file
            run.writeCal;
                        
            % Write input text file
            fid=fopen(fullfile(run.folder, 'input.txt'),'wt');
            fprintf(fid,' \n 0 \n');
            fclose(fid);
            
        end
                
        function [] = RunPreProcessor(run, varargin)
            % Runs the Nemoh preProcessor
            % By default, Matlab waits for Nemoh to finish running before
            % it continues.
            % 'Background' is an optional argument that runs Nemoh in the
            % background so that Matlab is free to use.
            
            opts = checkOptions({'Background'}, varargin);
            
            if (opts)
                system(['cd ' run.folder ' && ' fullfile(run.exePath, 'preProcessor'), '&']);
            else
                system(['cd ' run.folder ' && '  fullfile(run.exePath, 'preProcessor')]);
            end
        end
        
        function [] = RunSolver(run)
            % Runs the Nemoh Solver 
            % By default, Matlab waits for Nemoh to finish running before
            % it continues.
            % 'Background' is an optional argument that runs Nemoh in the
            % background so that Matlab is free to use.
            
            opts = checkOptions({'Background'}, varargin);
            
            if (opts)
                system(['cd ' run.folder ' && ' fullfile(run.exePath, 'Solver.exe'), '&']);
            else
                system(['cd ' run.folder ' && '  fullfile(run.exePath, 'Solver.exe')]);
            end
        end
        
        function [] = RunPostProcessor(run)
            % Runs the Nemoh postProcessor
            % By default, Matlab waits for Nemoh to finish running before
            % it continues.
            % 'Background' is an optional argument that runs Nemoh in the
            % background so that Matlab is free to use.
            
            opts = checkOptions({'Background'}, varargin);
            
            if (opts)
                system(['cd ' run.folder ' && ' fullfile(run.exePath, 'postProcessor.exe'), '&']);
            else
                system(['cd ' run.folder ' && '  fullfile(run.exePath, 'postProcessor.exe')]);
            end
        end
        
        function [] = Run(run, varargin)
            % Runs Nemoh (preProcessor, Solver, postProcessor) with system
            % command. 
            % By default, Matlab waits for Nemoh to finish running
            % before it continues. 'Background' is an optional argument
            % that runs Nemoh in the background so that Matlab is free to
            % use. 
            % By default, RunNemoh creates a batch file to run all of
            % preProcessor, Solver, postProcessor. 'NoBatch' is an optional
            % argument that does not create a batch file to run Nemoh.
            
            opts = checkOptions({{'Background'}, {'NoBatch'}, {'NoMesh'}}, varargin);
            back = opts(1);
            noBat = opts(2);
            noMesh = opts(3);
            
            if (~run.computeHS)
                noMesh = true;
            end
                        
            if (~noBat || back)
                filename = fullfile(run.folder,'nem_run.bat');
                fileID = fopen(filename, 'wt');
                if (~noMesh)
                    fprintf(fileID, '%s\\Mesh.exe >Mesh\\Mesh.log\n', run.exePath);
                end
                fprintf(fileID, '%s\\preProcessor.exe\n', run.exePath);
                fprintf(fileID, '%s\\Solver.exe\n', run.exePath);
                fprintf(fileID, '%s\\postProcessor.exe\n', run.exePath);
                fclose(fileID);
            end

            if (back)
                runStr = ['cd ' run.folder ' && nem_run.bat &'];
            else
                if (noBat)
                    runStr = ['cd ' run.folder];
                    if (~noMesh)
                        runStr = [runStr ' && ' run.exePath '/Mesh.exe >Mesh\Mesh.log '];
                    end
                    runStr = [runStr ' && ' run.exePath '/preProcessor.exe && ' run.exePath '/Solver.exe && ' run.exePath '/postProcessor.exe'];
                else
                    runStr = ['cd ' run.folder ' && nem_run.bat'];
                end
            end
            
            system(runStr);
        end
         
        function [Mass, Inertia, KH, XB, YB, ZB] = MakeAxisMesh(run, r, z, n, zcg, ntheta, Npan)
            
            %run.zcg = zcg;
            
            theta=[0.:pi/(ntheta-1):pi];
            nx=0;
            
            % Calcul des sommets du maillage
            for j=1:ntheta
                for i=1:n    
                    nx=nx+1;
                    x(nx)=r(i)*cos(theta(j));
                    y(nx)=r(i)*sin(theta(j));
                    z(nx)=z(i);
                end;
            end;
            
            % Calcul des facettes
            nf=0;
            for i=1:n-1
                for j=1:ntheta-1
                    nf=nf+1;
                    NN(1,nf)=i+n*(j-1);
                    NN(2,nf)=i+1+n*(j-1);
                    NN(3,nf)=i+1+n*j;
                    NN(4,nf)=i+n*j;
                end;
            end;
            
            % Affichage de la description du maillage
            nftri=0;
            for i=1:nf
                nftri=nftri+1;
                tri(nftri,:)=[NN(1,i) NN(2,i) NN(3,i)];
                nftri=nftri+1;
                tri(nftri,:)=[NN(1,i) NN(3,i) NN(4,i)];
            end;
            
            figure;
            trimesh(tri,x,y,z,[zeros(nx,1)]);
            axis equal;
            title('Characteristics of the discretisation');
            fprintf('\n --> Number of nodes             : %g',nx);
            fprintf('\n --> Number of panels (max 2000) : %g \n',nf);

            fid = fopen(fullfile(run.folder, 'mesh.cal'),'w');
            fprintf(fid, 'axisym \n');
            fprintf(fid,'1 \n 0. 0. \n ');
            
            fprintf(fid,'%f %f %f \n',[0. 0. zcg]);
            fprintf(fid,'%g \n 2 \n 0. \n 1.\n', Npan);
            
            fclose(fid);
            
            fid=fopen([run.folder '/Mesh/axisym'],'wt');
            fprintf(fid,'%g \n',nx);
            fprintf(fid,'%g \n',nf);
            for i=1:nx
                fprintf(fid,'%E %E %E \n',[x(i) y(i) z(i)]);
            end;
            
            for i=1:nf
                fprintf(fid,'%g %g %g %g \n',NN(:,i)');
            end;
            fclose(fid);
            
            Mass = 0;
            Inertia = 0; 
            KH = 0;
            XB = 0;
            YB = 0;
            ZB = 0;
            [Mass, Inertia, KH, XB, YB, ZB] = run.RunMesh;
        end
        
        function [Mass, Inertia, KH, XB, YB, ZB] = RunMesh(run)
            % Runs the Nemoh Mesh function
            oldFolder = cd(run.folder);
            system([fullfile(run.exePath,'Mesh.exe') '>' fullfile('Mesh','Mesh.log')]);
            cd(oldFolder);
            
            fid=fopen(fullfile(run.folder,'Mesh','axisym.tec'),'r');
            ligne=fscanf(fid,'%s',2);
            run.nx=fscanf(fid,'%g',1);
            ligne=fscanf(fid,'%s',2);
            run.nf=fscanf(fid,'%g',1);
            ligne=fgetl(fid);
            
            fprintf('\n Characteristics of the mesh for Nemoh \n');
            fprintf('\n --> Number of nodes : %g',run.nx);
            fprintf('\n --> Number of panels : %g\n \n',run.nf);
            
            for i=1:run.nx
                ligne=fscanf(fid,'%f',6);
                x(i)=ligne(1);
                y(i)=ligne(2);
                z(i)=ligne(3);
            end;
            
            for i=1:run.nf
                ligne=fscanf(fid,'%g',4);
                NN(1,i)=ligne(1);
                NN(2,i)=ligne(2);
                NN(3,i)=ligne(3);
                NN(4,i)=ligne(4);
            end;
            
            nftri=0;
            for i=1:run.nf
                nftri=nftri+1;
                tri(nftri,:)=[NN(1,i) NN(2,i) NN(3,i)];
                nftri=nftri+1;
                tri(nftri,:)=[NN(1,i) NN(3,i) NN(4,i)];
            end;

            ligne=fgetl(fid);
            ligne=fgetl(fid);
            for i=1:run.nf    
                ligne=fscanf(fid,'%g %g',6);
                xu(i)=ligne(1);
                yv(i)=ligne(2);
                zw(i)=ligne(3);
                u(i)=ligne(4);
                v(i)=ligne(5);
                w(i)=ligne(6);
            end;
            
            fclose(fid);
            figure;
            trimesh(tri,x,y,z);
            axis equal;
            hold on;
            quiver3(xu,yv,zw,u,v,w);
            title('Mesh for Nemoh');
            
            clear KH;
            KH = zeros(6,6);
            fid=fopen(fullfile(run.folder,'Mesh','KH.dat'),'r');
            
            for i=1:6   
                ligne=fscanf(fid,'%g %g',6);
                KH(i,:)=ligne;
            end;
            
            fclose(fid);
            clear XB YB ZB Mass WPA Inertia
            
            Inertia=zeros(6,6);
            fid=fopen(fullfile(run.folder,'Mesh','Hydrostatics.dat'),'r');
            ligne=fscanf(fid,'%s',2);
            XB=fscanf(fid,'%f',1);
            ligne=fgetl(fid);
            ligne=fscanf(fid,'%s',2);
            YB=fscanf(fid,'%f',1);
            ligne=fgetl(fid);
            ligne=fscanf(fid,'%s',2);
            ZB=fscanf(fid,'%f',1);
            ligne=fgetl(fid);
            ligne=fscanf(fid,'%s',2);
            Mass=fscanf(fid,'%f',1)*run.rho;
            ligne=fgetl(fid);
            ligne=fscanf(fid,'%s',2);
            WPA=fscanf(fid,'%f',1);
            status=fclose(fid);
            clear ligne
            
            fid=fopen(fullfile(run.folder,'Mesh','Inertia_hull.dat'),'r');
            
            for i=1:3
                ligne=fscanf(fid,'%g %g',3);
                Inertia(i+3,4:6)=ligne;
            end;
            
            Inertia(1,1)=Mass;
            Inertia(2,2)=Mass;
            Inertia(3,3)=Mass;
        end
    end
    
    methods (Access = protected)
        function [] = makeFolderAndSub(run, fold)
            
            if (~exist(fold, 'dir'))
                system(['mkdir ', fold]);
            end
            
            if (~exist(fullfile(fold,'Mesh'), 'dir'))
                system(['mkdir ', fullfile(fold, 'Mesh')]);
            end
            
            if (~exist(fullfile(fold, 'Results'), 'dir'))
                system(['mkdir ', fullfile(fold, 'Results')]);
            end
            
            % Write ID file
            fid = fopen(fullfile(run.folder, 'ID.dat'),'wt');
            fprintf(fid, '1\n.');
            fclose(fid);
        end
    end
    
    methods (Access = private)
        function [] = writeCal(run)
            if (run.nOmega < 1)
                error('The number of frequencies must be greater than 0');
            end
            
            if (run.nBeta < 1)
                error('The number of directions must be greater than 0');
            end
            
            if (isempty(run.omegaLim))
                error('It does not appear that the frequency limits have been set.');
            end
            
            if (isempty(run.betaLim))
                error('It does not appear that the direction limits have been set.');
            end
            
            if (isempty(run.h))
                error('Depth has not been set');
            end
            
            if (isempty(run.floatBods))
                error('It does not appear that the Floating bodies have been set');
            end
            
            nbod = length(run.floatBods);
            
            run.geoFiles = cell(1, nbod);
            for n = 1:nbod
                if (isempty(run.floatBods(n).Handle))
                    name = ['body'  num2str(n)];
                else
                    name = run.floatBods(n).Handle;
                end
                geoFile = run.writeMeshFile(run.floatBods(n), false, name);
                run.geoFiles{n} = geoFile;
            end
            
            % write Cal file
            fid = fopen(fullfile(run.folder, 'Nemoh.cal'),'wt');
                      
            fprintf(fid,'--- Environment ------------------------------------------------------------------------------------------------------------------ \n');
            fprintf(fid,'%7.2f                                              ! RHO               ! KG/M**3       ! Fluid specific volume \n', run.rho);
            fprintf(fid,'%7.2f                                              ! G                 ! M/S**2        ! Gravity \n', IWaves.G);
            fprintf(fid,'%7.2f                                              ! DEPTH             ! M             ! Water depth\n', run.h);
            fprintf(fid,'0.     0.                                            ! XEFF YEFF         ! M             ! Wave measurement point\n');
            fprintf(fid,'--- Description of floating bodies -----------------------------------------------------------------------------------------------\n');
            fprintf(fid,'%4i                                                 ! Number of bodies\n', nbod);
            
            for n = 1:nbod
                run.writeBodyCal(fid, n, run.floatBods(n), run.geoFiles{n});
            end
            betad1 = round(1e4*180/pi*run.betaLim(1))/1e4;
            betad2 = round(1e4*180/pi*run.betaLim(2))/1e4;
            fprintf(fid,'--- Load cases to be solved -------------------------------------------------------------------------------------------------------\n');
            fprintf(fid,'%4i     %7.4f   %7.4f                           ! Number of wave frequencies, Min, and Max (rad/s)\n', run.nOmega, run.omegaLim(1), run.omegaLim(2));
            fprintf(fid,'%4i     %7.2f   %7.2f                           ! Number of wave directions, Min and Max (degrees)\n', run.nBeta, betad1, betad2);
            fprintf(fid,'--- Post processing ---------------------------------------------------------------------------------------------------------------\n');
            fprintf(fid,'1  0.1     10.                                       ! IRF                     ! IRF calculation (0 for no calc), time step and duration\n');
            fprintf(fid,'0                                                    ! Show pressure\n');
            fprintf(fid,'0	0.	180.                                 ! Kochin function         ! Number of directions of calculation (0 for no calc), Min and Max (degrees)\n');
            
            if (~isempty(run.fieldArray))
                nps = run.fieldArray.NumberPoints;
                lens = run.fieldArray.Lengths;
            else
                nps = [0 0];
                lens = [0 0 0];
            end
            fprintf(fid,'%5i    %5i  %7.2f   %7.2f                    ! Free surface elevation  ! Number of points in x direction (0 for no calc) and y direction and dimensions of domain in x and y direction\n', nps(1), nps(2), lens(1), lens(2));	
            if (run.useCylSurf)
                if (~isempty(run.cylArray))
                    r = run.cylArray.Radius;
                    nth = run.cylArray.Ntheta;
                    nz = run.cylArray.Nz;
                else
                    r = 0;
                    nth = 1;
                    nz = 1;
                end
                fprintf(fid,'%7.3f  %5i  %5i                                ! Cylindrical surface     ! Radius (0 for no calc), Number of points in theta direction, Number of points in z direction\n', r, nth, nz);	
            end
            fprintf(fid,'--- End of file -------------------------------------------------------------------------------------------------------------------\n');
            
            fclose(fid);
        end
        
        function [] = writeBodyCal(run, fid, n, floatBod, geoFile)
            pos = [floatBod.XYpos floatBod.Zpos] + floatBod.CenterRot;
            modes = floatBod.Modes;
            vector = modes.Vector;
            
            fprintf(fid,'--- Body %i -----------------------------------------------------------------------------------------------------------------------\n', n);
            fprintf(fid,'%s\\Mesh\\%s		! Name of mesh file\n', run.folder, geoFile);
            fprintf(fid,'%6i    %6i                                     ! Number of points and number of panels 	\n', run.nNode, run.nPan);
            fprintf(fid,'%4i                                                 ! Number of degrees of freedom\n', modes.DoF);
            
            if (vector(1))
                fprintf(fid,'1 1. 0. 0. %7.3f %7.3f %7.3f                   ! Surge\n', pos);
            end
            if (vector(2))
                fprintf(fid,'1 0. 1. 0. %7.3f %7.3f %7.3f                   ! Sway\n', pos);
            end
            if (vector(3))
                fprintf(fid,'1 0. 0. 1. %7.3f %7.3f %7.3f                   ! Heave\n', pos);
            end
            if (vector(4))
                fprintf(fid,'2 1. 0. 0. %7.3f %7.3f %7.3f                   ! Roll about a point\n', pos);
            end
            if (vector(5))
                fprintf(fid,'2 0. 1. 0. %7.3f %7.3f %7.3f                   ! Pitch about a point\n', pos);
            end
            if (vector(6))
                fprintf(fid,'2 0. 0. 1. %7.3f %7.3f %7.3f                   ! Yaw about a point\n', pos);
            end
            
            fprintf(fid,'%4i                                                 ! Number of resulting generalised forces\n', modes.DoF);
            if (vector(1))
                fprintf(fid,'1 1. 0. 0. %7.3f %7.3f %7.3f                   ! Force in x direction\n', pos);
            end
            if (vector(2))
                fprintf(fid,'1 0. 1. 0. %7.3f %7.3f %7.3f                   ! Force in y direction\n', pos);
            end
            if (vector(3))
                fprintf(fid,'1 0. 0. 1. %7.3f %7.3f %7.3f                   ! Force in z direction\n', pos);
            end
            if (vector(4))
                fprintf(fid,'2 1. 0. 0. %7.3f %7.3f %7.3f                   ! Moment force in x direction about a point\n', pos);
            end
            if (vector(5))
                fprintf(fid,'2 0. 1. 0. %7.3f %7.3f %7.3f                   ! Moment force in y direction about a point\n', pos);
            end
            if (vector(6))
                fprintf(fid,'2 0. 0. 1. %7.3f %7.3f %7.3f                   ! Moment force in z direction about a point\n', pos);
            end
            
            fprintf(fid,'0                                                    ! Number of lines of additional information \n');
            
        end
        
        function [] = writeMeshCalAndFile(run, floatBod)
            
            if (isempty(floatBod.Handle))
                name = 'body1';
            else
                name = floatBod.Handle;
            end
            geoFile = run.writeMeshFile(floatBod, true, name);
            
            fid = fopen(fullfile(run.folder, 'mesh.cal'),'wt');
            fprintf(fid, '%s \n', geoFile);
            fprintf(fid, '1 \n 0. 0. \n ');
            
            fprintf(fid, '%f %f %f \n', floatBod.Cg);
            geo = PanelGeo(floatBod.PanelGeo);
            npan = sum((geo.IsWets + ~geo.IsInteriors) == 2);
            fprintf(fid,'%g \n 2 \n 0. \n 1.\n', npan);
            
            fclose(fid);
        end
        
        function [geoFile] = writeMeshFile(run, floatBod, forMesh, name)
            geo = PanelGeo(floatBod.PanelGeo);
            
            if (geo.Ysymmetry)
                error('Body Y-symmetry not allowed in Nemoh');
            end

            pos = [floatBod.XYpos floatBod.Zpos];
            geo.Translate(pos);
            
            geoPans = geo.Panels;
            npanTot = geo.Count;
            % only include panels that are wet and are not the interior
            % surface
            npan = sum((geo.IsWets + ~geo.IsInteriors) == 2);
            
            % write it the simple way first (multiples of same nodes)
            nodes = zeros(npan*4, 3);
            pans = zeros(npan, 4);
            
            m = 0;
            for nn = 1:npanTot
                pan = geoPans(nn);
                if (pan.IsWet && ~pan.IsInterior)
                    m = m + 1;
                    verts = pan.Vertices;
                    nstart = (m - 1)*4;

                    for n = 1:4
                        inode = nstart + n;
                        nodes(inode,:) = verts(n,:);
                        pans(m, n) = inode;
                    end
                end
            end
            
            geoFile = [name '.dat'];
            
            fid = fopen(fullfile(run.folder, 'Mesh', geoFile),'wt');
            
            if (forMesh)
                fprintf(fid,' %i \n', npan*4);
                fprintf(fid,' %i \n', npan);
            else
                fprintf(fid,' 2 %i \n', geo.Xsymmetry);
            end
            
            for n = 1:npan*4
                if (forMesh)
                    fprintf(fid, '%8.4f %8.4f %8.4f \n', nodes(n,:));
                else
                    fprintf(fid, '%i %8.4f %8.4f %8.4f \n', n, nodes(n,:));
                end
            end
            if (~forMesh)
                fprintf(fid, '0 0.0 0.0 0.0 \n');
            end
            for n = 1:npan
                fprintf(fid, '%i %i %i %i \n', pans(n,:));
            end
            
            fclose(fid);
            
            run.nNode = size(nodes, 1);
            run.nPan = size(pans, 1);
        end
    end
end