% This test sets up and runs Wamit on a 6 dof cylinder. 
% The mean drift forces are computed using a control surface, as in 
% example_Wamit_driftForces_1. 
% The difference is that this script computes the body motions
% outside of WAMIT and feeds the RAOs back in for computation of the drift
% forces. This may be of limited use for a single body, but for devices with
% constraints not handled by WAMIT, the drift forces can still be computed
% accurately using WAMIT.

clear; close all; clc;

%% Set up the run

% every run needs a name. All the Wamit input and output files are named 
% based on this name
run_name = 'wam_1b_6dof_4';         

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
cyl.Modes = ModesOfMotion([1 0 1 0 1 0]);   % This sets up the modes of 
                                            % motion for the cylinder.
                                            
% Create a WamitRunCondition object - this builds the Wamit run
wam_run = WamitRunCondition(folder, run_name);  

wam_run.Rho = rho;      % set the fluid density (kg/m^3)
wam_run.T = 4:0.1:12;   % set the wave period(s)(in s)
wam_run.Beta = 45*pi/180;       % set the incident wave direction(s) (in radians)
wam_run.H = Inf;        % set the water depth (in m) - can be a positive 
                        % value or Inf
wam_run.CompDrift = true; % Turn drift force computation on/off
wam_run.NCPU = 6; % Set no. of processors to speed things up a bit
wam_run.RAMGBmax = 16; % Set amount of RAM available (though this example 
                      % is unlikely to exceed even a smaller limit).

if wam_run.CompDrift
    wam_run.DriftOption = [1 2 3]; % Choose methods for computing drift forces (can choose more than one) 
                            % 1 - using control surface(s) [IOPTN(7) in .FRC file] 
                            % 2 - momentum conservation [IOPTN(8) in .FRC file]
                            % 3 - pressure integration [IOPTN(9) in .FRC file]
    if ismember(1,wam_run.DriftOption)
        wam_run.AutoCSF = true; % Decide whether to use user-input control surface (false) 
                         % or to use the automatic creation feature of WAMIT (true).
        wam_run.BoxCSF = false; % Use cylinder or quadrilateral-based box shapes control surface (true - box, false - cylinder)
    end
end

wam_run.IReadRAO = 1;               % 0 - use .4 file RAOs, 1 - use user-input RAOs.

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

% Compute RAOs externally to WAMIT
% 1) Get the frequency domain forces computed by WAMIT.
% i.e. the added mass, hydrodynamic damping, excitation force, and
% hydrostatic stiffness matrix
freqDomForces = wam_result.FreqDomForces;

% 2) Create the constraint matrix to constrain the bodies together.
P_6dof = ConstraintMatComp.FixedBodies([0 0 0]); % This should simply give the identity matrix in the case of a single body
P = P_6dof(logical(cyl.Modes.Vector),...
logical(cyl.Modes.Vector)); % Remove rows and columns for DoFs not studied.

% 3) Create the Frequency Domain computation - type of IEnergyComp
freqDomComp = FreqDomComp(freqDomForces, cyl, 'Constrained', P);

% 4) Add desired changes to damping and stiffness matrices
Dpto_6dof = zeros(6,6); Dpto_6dof(3,3) = 100000; 
Dpto = Dpto_6dof(logical(cyl.Modes.Vector),...
    logical(cyl.Modes.Vector));
freqDomComp.SetDpto(Dpto);

% 5) Derive motion RAOs
motionRAOs = freqDomComp.Motions;

%% Write RAOs to .rao file for input to WAMIT and alter .CFG file
wam_run.MotionRAOs = motionRAOs;
wam_run.WriteRunIReadRAO;

% Rerun WAMIT to obtain drift forces
wam_run.Run;

%% Read results and plot drift forces
wam_result = WamitResult(wam_run);  % create the updated WamitRunCondition
wam_result.ReadResult; % Read the updated WAMIT result

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
 
