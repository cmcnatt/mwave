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
classdef NemohRunCondition < handle
    % Creates input files for a Nemoh run
    
    properties (Access = private)
        nemohPath;
        folder;
        zcg;
        nx;
        nf;
        rho;
        h;
        omegaLim;
        nOmega;
        betaLim;
        nBeta;
    end

    properties (Dependent)
        Rho;                % Fluid density
        H;                  % Water depth
        OmegaLimits;        % Lower and upper limit of the radial frequencies
        OmegaCount;         % Number of radial frequencies
        BetaLimits;         % Lower and upper limit of the wave direction
        BetaCount;          % Number of wave directions
        T;                  % The wave periods (can't set)
        Beta;               % The wave directions (can't set)
        NemohPath;          % Path location of Nemoh executables
        Folder;             % Folder location of run (string)
    end
    
    methods
        
        function [run] = NemohRunCondition(folder)
            run.nemohPath = [myMatPath '\..\nemoh'];
            run.rho = 1025;
            if (nargin == 0)
                run.folder = ' ';
            else 
                run.folder = folder;
                NemohRunCondition.makeFolderAndSub(folder);
            end
        end
        
        function [fol] = get.Folder(run)
            % Get the folder location where the input files will go
            fol = run.folder;
        end
        function [run] = set.Folder(run, fold)
            % Set the folder location where the input files will go
            if (ischar(fold))
                run.folder = fold;
                NemohRunCondition.makeFolderAndSub(fold);
            else
                error('Folder must be a string');
            end
        end
        
        function [np] = get.NemohPath(run)
            % Path location of nemoh executables
            np = run.nemohPath;
        end
        function [run] = set.NemohPath(run, np)
            % Path location of nemoh executables
            if (ischar(np))
                run.nemohPath = np;
            else
                error('NemohPath must be a string');
            end
        end
        
        function [rh] = get.Rho(run)
            % Get the water density
            rh = run.rho;
        end
        function [run] = set.Rho(run, rh)
            % Set the water density
            if (rh > 0)
                run.rho = rh;
            else
                error('Rho must be a positive number');
            end
        end
        
        function [h_] = get.H(run)
            % Get the water depeth
            h_ = run.h;
        end
        function [run] = set.H(run, h_)
            % Set the water depth
            if (h_ > 0)
                run.h = h_;
            else
                error('The depth must be a positive number');
            end
        end
        
        function [ol] = get.OmegaLimits(run)
            % Get the lower and upper freguency limits
            ol = run.omegaLim;
        end
        function [run] = set.OmegaLimits(run, ol)
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
        
        function [ol] = get.BetaLimits(run)
            % Get the lower and upper incident direction limits in degrees
            ol = run.omegaLim;
        end
        function [run] = set.BetaLimits(run, bl)
            % Get the lower and upper incidnet direction limits in degrees
            
            if (length(bl) ~= 2)
                error('The limits must be a vector of length 2');
            end
            
            if (bl(2) < bl(1))
                error('The first element is the lower limit, which must be less than the upper limit');
            end
            
            if any(bl <= 0)
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
        
        function [Mass, Inertia, KH, XB, YB, ZB] = MakeAxisMesh(run, r, z, n, zcg, ntheta, Npan)
            
            run.zcg = zcg;
            
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

            fid = fopen([run.folder '\mesh.cal'],'w');
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
            
            [Mass, Inertia, KH, XB, YB, ZB] = run.RunMesh;
        end
        
        function [] = WriteNemohCal(run)
            fid = fopen([run.folder, '\Nemoh.cal'],'wt');
            
            fprintf(fid,'--- Environment ------------------------------------------------------------------------------------------------------------------ \n');
            fprintf(fid,'%f				! RHO 			! KG/M**3 	! Fluid specific volume \n', run.rho);
            fprintf(fid,'9.81				! G			! M/S**2	! Gravity \n');
            fprintf(fid,'0.                 ! DEPTH			! M		! Water depth\n');
            fprintf(fid,'0.	0.              ! XEFF YEFF		! M		! Wave measurement point\n');
            fprintf(fid,'--- Description of floating bodies -----------------------------------------------------------------------------------------------\n');
            fprintf(fid,'1				! Number of bodies\n');
            fprintf(fid,'--- Body 1 -----------------------------------------------------------------------------------------------------------------------\n');
            fprintf(fid,'%s\\Mesh\\axisym.dat		! Name of mesh file\n', run.folder);
            fprintf(fid,'%g %g			! Number of points and number of panels 	\n', run.nx, run.nf);
            fprintf(fid,'6				! Number of degrees of freedom\n');
            fprintf(fid,'1 1. 0.	0. 0. 0. 0.		! Surge\n');
            fprintf(fid,'1 0. 1.	0. 0. 0. 0.		! Sway\n');
            fprintf(fid,'1 0. 0. 1. 0. 0. 0.		! Heave\n');
            fprintf(fid,'2 1. 0. 0. 0. 0. %f		! Roll about a point\n',run.zcg);
            fprintf(fid,'2 0. 1. 0. 0. 0. %f		! Pitch about a point\n',run.zcg);
            fprintf(fid,'2 0. 0. 1. 0. 0. %f		! Yaw about a point\n',run.zcg);
            fprintf(fid,'6				! Number of resulting generalised forces\n');
            fprintf(fid,'1 1. 0.	0. 0. 0. 0.		! Force in x direction\n');
            fprintf(fid,'1 0. 1.	0. 0. 0. 0.		! Force in y direction\n');
            fprintf(fid,'1 0. 0. 1. 0. 0. 0.		! Force in z direction\n');
            fprintf(fid,'2 1. 0. 0. 0. 0. %f		! Moment force in x direction about a point\n',run.zcg);
            fprintf(fid,'2 0. 1. 0. 0. 0. %f		! Moment force in y direction about a point\n',run.zcg);
            fprintf(fid,'2 0. 0. 1. 0. 0. %f		! Moment force in z direction about a point\n',run.zcg);
            fprintf(fid,'0				! Number of lines of additional information \n');
            fprintf(fid,'--- Load cases to be solved -------------------------------------------------------------------------------------------------------\n');
            fprintf(fid,'1	0.8	0.8		! Number of wave frequencies, Min, and Max (rad/s)\n');
            fprintf(fid,'1	0.	0.		! Number of wave directions, Min and Max (degrees)\n');
            fprintf(fid,'--- Post processing ---------------------------------------------------------------------------------------------------------------\n');
            fprintf(fid,'1	0.1	10.		! IRF 				! IRF calculation (0 for no calculation), time step and duration\n');
            fprintf(fid,'0				! Show pressure\n');
            fprintf(fid,'0	0.	180.		! Kochin function 		! Number of directions of calculation (0 for no calculations), Min and Max (degrees)\n');
            fprintf(fid,'0	50	400.	400.	! Free surface elevation 	! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction\n');	
            fprintf(fid,'---');
            
            fclose(fid);
        end
        
        function [Mass, Inertia, KH, XB, YB, ZB] = RunMesh(run)
            % Runs the Nemoh Mesh function
            oldFolder = cd(run.folder);
            system([run.nemohPath '\Mesh.exe >Mesh\Mesh.log']);
            cd(oldFolder);
            
            fid=fopen([run.folder,'\Mesh\axisym.tec'],'r');
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
            fid=fopen([run.folder,'\Mesh\KH.dat'],'r');
            
            for i=1:6   
                ligne=fscanf(fid,'%g %g',6);
                KH(i,:)=ligne;
            end;
            
            fclose(fid);
            clear XB YB ZB Mass WPA Inertia
            
            Inertia=zeros(6,6);
            fid=fopen([run.folder,'\Mesh\Hydrostatics.dat'],'r');
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
            
            fid=fopen([run.folder,'\Mesh\Inertia_hull.dat'],'r');
            
            for i=1:3
                ligne=fscanf(fid,'%g %g',3);
                Inertia(i+3,4:6)=ligne;
            end;
            
            Inertia(1,1)=Mass;
            Inertia(2,2)=Mass;
            Inertia(3,3)=Mass;
        end
        
        function [] = RunPreProcessor(run)
            % Runs the Nemoh preProcessor
            % By default, Matlab waits for Nemoh to finish running before
            % it continues.
            % 'Background' is an optional argument that runs Nemoh in the
            % background so that Matlab is free to use.
            
            opts = checkOptions({'Background'}, varargin);
            
            if (opts)
                system(['cd ' run.folder ' && ' run.nemohPath '\preProcessor.exe &']);
            else
                system(['cd ' run.folder ' && '  run.nemohPath '\preProcessor.exe']);
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
                system(['cd ' run.folder ' && ' run.nemohPath '\Solver.exe &']);
            else
                system(['cd ' run.folder ' && '  run.nemohPath '\Solver.exe']);
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
                system(['cd ' run.folder ' && ' run.nemohPath '\postProcessor.exe &']);
            else
                system(['cd ' run.folder ' && '  run.nemohPath '\postProcessor.exe']);
            end
        end
        
        function [] = RunNemoh(run, w, dir, depth, varargin)
            % Runs Nemoh (preProcessor, Solver, postProcessor) with system
            % command. 
            % By default, Matlab waits for Nemoh to finish running
            % before it continues. 'Background' is an optional argument
            % that runs Nemoh in the background so that Matlab is free to
            % use. 
            % By default, RunNemoh creates a batch file to run all of
            % preProcessor, Solver, postProcessor. 'NoBatch' is an optional
            % argument that does not create a batch file to run Nemoh.
            
            opts = checkOptions({{'Background'}, {'NoBatch'}}, varargin);
            back = opts(1);
            noBat = opts(2);
            
            fid=fopen([run.folder,'\Nemoh.cal'],'r');
            for i=1:6
                ligne = fgetl(fid);
            end
            nBodies = fscanf(fid,'%g',1);
            fclose(fid);
            
            fid=fopen([run.folder,'\Nemoh.cal'],'r');
            n=1;
            clear textline;
            textline={};
            while (~feof(fid))
                textline(n)={fgetl(fid)};
                if (n == 4) 
                    textline(n)={sprintf('%f                 ! DEPTH			! M		! Water depth',depth)};
                end
                if ((mod(n,18) == 9) && ((n-9)/18 <= nBodies))
                    temp=cell2mat(textline(n));
                    temp2=[];
                    ntemp=length(temp);
                    k=1;
                    for i=1:ntemp
                        if (temp(i) == '\')
                            temp2=[temp2,temp(k:i),'\'];
                            k=i+1;
                        end;            
                    end
                    temp2=[temp2,temp(k:ntemp)];
                    textline(n)={temp2};
                    cell2mat(textline(n));
                end
                if (n == 9+18*nBodies)
                    textline(n)={sprintf('%g %f %f       ! Number of wave frequencies, Min, and Max (rad/s)',length(w),w(1),w(length(w)))};
                end
                if (n == 10+18*nBodies)
                    textline(n)={sprintf('%g %f %f		! Number of wave directions, Min and Max (degrees)',1,dir,dir)};
                end
                n=n+1;
            end
            fclose(fid);
            
            fid = fopen([run.folder, '\Nemoh.cal'], 'wt'); 
            for i=1:n-1
                fprintf(fid, [cell2mat(textline(i)),'\n']);
            end
            fclose(fid);
            fid=fopen([run.folder, '/input.txt'],'wt');
            fprintf(fid,' \n 0 \n');
            fclose(fid);
            
            if (~noBat || back)
                filename = [run.folder '\nem_run.bat'];
                fileID = fopen(filename, 'wt');
                fprintf(fileID, '%s\\preProcessor.exe\n', run.nemohPath);
                fprintf(fileID, '%s\\Solver.exe\n', run.nemohPath);
                fprintf(fileID, '%s\\postProcessor.exe\n', run.nemohPath);
                fclose(fileID);
            end
                            
            % Calcul des coefficients hydrodynamiques
            if (back)
                system(['cd ' run.folder ' && nem_run.bat &']);
            else
                if (noBat)
                    system(['cd ' run.folder ' && ' run.nemohPath '\preProcessor.exe && ' run.nemohPath '\Solver.exe && ' run.nemohPath '\postProcessor.exe'])
                else
                    system(['cd ' run.folder ' && nem_run.bat']);
                end
            end
        end
        
         function [A, B, Fe] = ReadNemoh(run)
             % Lecture des resultats CA CM Fe
             clear Periode A B Famp Fphi Fe;
             fid = fopen([run.folder,'\Nemoh.cal'],'r');
             for n = 1:6
                 fgetl(fid);
             end
             nBodies=fscanf(fid,'%g',1);
             for n = 1:(2 + 18*nBodies)
                 fgetl(fid);
             end
             nw=fscanf(fid,'%g',1);
             fclose(fid);
             
             fid = fopen([run.folder,'\Results\ExcitationForce.tec'],'r');
             fgetl(fid);
             for c=1:6*nBodies
                 fgetl(fid);
             end;
             
             fgetl(fid);
             for k=1:nw
                 ligne = fscanf(fid,'%f',1+12*nBodies);
                 w(n) = ligne(1);
                 for j=1:6*nBodies
                     Famp(k,j) = ligne(2*j);
                     Fphi(k,j) = ligne(2*j+1);
                 end;
             end;
             
             fclose(fid);
             
             fid = fopen([run.folder,'\Results\RadiationCoefficients.tec'],'r');
             fgetl(fid);
             for n=1:6*nBodies
                fgetl(fid);
             end;
             
             for n=1:nBodies*6
                 fgetl(fid);
                 for k=1:nw
                     ligne=fscanf(fid,'%f',1+12*nBodies);
                     for j=1:6*nBodies
                         A(n,j,k)=ligne(2*j);
                         B(n,j,k)=ligne(2*j+1);
                     end;
                     fgetl(fid);
                 end;
             end;
             fclose(fid);
             
             % Expression des efforts d excitation de houle sous la forme Famp*exp(i*Fphi)
             Fe = Famp.*(cos(Fphi)+1i*sin(Fphi));
             
         end
    end
    
    methods (Static, Access = private)
        function [] = makeFolderAndSub(fold)
            
            if (~exist(fold, 'dir'))
                system(['mkdir ', fold]);
            end
            
            if (~exist([fold '\Mesh'], 'dir'))
                system(['mkdir ', fold, '\Mesh']);
            end
            
            if (~exist([fold '\Results'], 'dir'))
                system(['mkdir ', fold, '\Results']);
            end
            
            fid = fopen([fold '\ID.dat'],'wt');
            fprintf(fid, '1\n.');
            fclose(fid);
        end
    end
end