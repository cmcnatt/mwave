% This test sets up and runs Wamit on a 6 dof cylinder and the mean drift 
% forces are computed using a control surface in WAMIT. The
% cylinder is the same as example_Wamit_1bod_6dof_2

clear; close all; clc;

%% Set up the run

% every run needs a name. All the Wamit input and output files are named 
% based on this name
run_name = 'wam_1b_6dof_3';         

% This is the folder where all the Wamit input and output files go 
folder = [mwavePath '\Examples\bemRunFolder'];  
if 0 == exist(folder, 'dir')
    mkdir(folder);
end
delete([folder '\*']);
                                                    
% mwavePath is a function that returns the path to the to the folder where
% mwave is kept it locates this by finding a file called mwave.home     

rho = 1000;     % the fluid density (kg/m^2). Set in the Wamit run and used
                % to compute body properties

% A run needs a floating body. In this case it is a cylinder. Here, we're
% going to build the floating body with an mwave class.

% The cylinder parameters are:
diameter = 10;
draft = 10;
height = 15;            % 10 m below the free surface, 5 m above

% These parameters described the number of panels that will be used to
% create the panel mesh (.gdf file)
Ntheta = 36;            % Number in the circular direction
Nr = 6;                 % Number in the radial direction
Nz = 18;                % Number in the vertical direction

cyl = FloatingCylinder(rho, diameter/2, height, draft, Ntheta, Nr, Nz);

% We can look at the Geometry..
figure;
plot(cyl.PanelGeo);
axis equal;
box on;
grid on;
set(gca, 'view', [-10 22]);
title({'Cylinder Panelization', 'dia: 10 m, draft: 10 m, height: 15 m',});

% Like in  wam_1b_6dof_1, we still need to do these things, and the rest of
% the computation is the same as wam_1b_6dof_1
cyl.Handle = 'Cylinder';        % Set its name. This is more relavent with 
                                % multiple bodies                                    
cyl.Modes = ModesOfMotion([1 1 1 1 1 1]);   % This sets up the modes of 
                                            % motion for the cylinder.
                                            
% Create a WamitRunCondition object - this builds the Wamit run
wam_run = WamitRunCondition(folder, run_name);  

wam_run.Rho = rho;      % set the fluid density (kg/m^3)
wam_run.T = 4:0.1:12;   % set the wave period(s)(in s)
wam_run.Beta = 0;       % set the incident wave direction(s) (in radians)
wam_run.H = Inf;        % set the water depth (in m) - can be a positive 
                        % value or Inf
wam_run.CompDrift = true; % Turn drift force computation on/off
wam_run.NCPU = 6;       % Set no. of processors to speed things up a bit
wam_run.RAMGBmax = 16;  % Set amount of RAM available (though this example 
                        % is unlikely to exceed even a smaller limit).

if wam_run.CompDrift
    wam_run.DriftOption = [1 2 3];  % Choose methods for computing drift forces (can choose more than one) 
                                    % 1 - using control surface(s) [IOPTN(7) in .FRC file] 
                                    % 2 - momentum conservation [IOPTN(8) in .FRC file]
                                    % 3 - pressure integration [IOPTN(9) in .FRC file]
    if ismember(1, wam_run.DriftOption)
        wam_run.AutoCSF = true; % Decide whether to use user-input control surface (false) 
                                % or to use the automatic creation feature of WAMIT (true).
        wam_run.BoxCSF = true;  % Use cylinder or quadrilateral-based box shapes control surface (true - box, false - cylinder)
    end
end

wam_run.FloatingBodies = cyl;       % Now set the floating body we just 
                                    % created as the body used in Wamit run
wam_run.WriteRun;                   % this writes all the necessary Wamit 
                                    % input files
                                    
%% Run Wamit

wam_run.Run;                            % Runs Wamit on the run that was 
                                        % just created. 
                                        
%% Read results

% the WamitResult class reads in Wamit output files. 
wam_result = WamitResult(wam_run);  % create with the WamitRunCondition

wam_result.ReadResult;              % Read the Wamit output files

driftForces = wam_result.DriftForces; % Extract drift forces

% Plot drift forces for all 6 DoFs of the cylinder and for all methods
legendNames = {'Control Surface','Momentum','Pressure Integration'};
dofsTitle = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};

% Plot different methods together on same figures
figure
for i = 1:6
    subplot(3,2,i)
    hold on
    for j = 1:size(driftForces,2)
        plot(driftForces(j).T,abs(driftForces(j).meanDriftForces(:,i)))
    end
    xlabel('Wave period / s')
    ylabel('Force/Moment')
    title(dofsTitle{i})
    grid on
    box on
end
legend(legendNames)

%% Plot different methods separately on different figures

for j = 1:size(driftForces,2)
    figure
    for i = 1:6
        subplot(3,2,i)
        plot(driftForces(j).T,abs(driftForces(j).meanDriftForces(:,i)))
        xlabel('Wave period / s')
        ylabel('Force/Moment')
        title(dofsTitle{i})
        grid on
        box on
    end
    sgtitle(driftForces(j).methodName)
end
