function [comp] = floatingCylComp(varargin)
            
    run = false;
    if ~isempty(varargin)
        run = varargin{1};
    end

    saveFile = [mwavePath 'QuickComps\matFiles\floatingCyl'];

    if run || ~exist([saveFile '.mat'])
        
        run_name = 'wam_1b_6dof_2';         
        folder = [mwavePath '\Examples\BemRuns\' run_name];  
        rho = 1025;     
        
        diameter = 10;
        draft = 10;
        height = 15;            
        
        Ntheta = 48;            % Number in the circular direction
        Nr = 8;                 % Number in the radial direction
        Nz = 24;                % Number in the vertical direction

        cyl = FloatingCylinder(rho, diameter/2, height, draft, Ntheta, Nr, Nz);
        
        cyl.Handle = 'Cylinder';                                       
        cyl.Modes = ModesOfMotion([1 1 1 1 1 1]);   
        wam_run = WamitRunCondition(folder, run_name);  

        wam_run.Rho = rho;      % set the fluid density (kg/m^3)
        wam_run.T = 4:0.1:12;   % set the wave period(s)(in s)
        wam_run.Beta = 0;       % set the incident wave direction(s) (in radians)
        wam_run.H = Inf;        % set the water depth (in m) - can be a positive 
                                % value or Inf     

        wam_run.FloatingBodies = cyl;       
        wam_run.WriteRun;                   
        wam_run.Run;                           
        wam_result = WamitResult(wam_run);  
        wam_result.ReadResult;              

        forces = wam_result.FreqDomForces;

        comp = FreqDomComp(forces, cyl);

        dof = cyl.Modes.DoF;
        Dpto = zeros(dof, dof);     
        Dpto(3,3) = 10^5;           
        comp.SetDpto(Dpto);    

        save(saveFile, 'comp');
    else
        load(saveFile);
    end

end