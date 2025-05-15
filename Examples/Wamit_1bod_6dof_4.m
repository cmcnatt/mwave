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
    C. McNatt, A. Cotten
%}
% This test sets up and runs Wamit on a 6 dof cylinder. 
% The mean drift forces are computed using a control surface, as in 
% Wamit_1bod_6dof_3. 
% The difference is that this script computes the body motions
% outside of WAMIT and feeds the RAOs back in for computation of the drift
% forces. This may be of limited use for a single body, but for devices with
% constraints not handled by WAMIT, the drift forces can still be computed
% accurately using WAMIT.

%% Set up the run

% every run needs a name. All the Wamit input and output files are named 
% based on this name
run_name = 'wam_1b_6dof_4';         

% This is the folder where all the Wamit input and output files go. 
% It is good practice to name the folder the 
% same as the run name.
folder = [mwavePath '\Examples\BemRuns\' run_name];
if ~exist(folder)
    mkdir(folder)
end
                                                    
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
Ntheta = 48;            % Number in the circular direction
Nr = 8;                 % Number in the radial direction
Nz = 24;                % Number in the vertical direction

cyl = FloatingCylinder(rho, diameter/2, height, draft, Ntheta, Nr, Nz);

% The FloatingCylinder class computes all of the mass properties that we
% had to compute by hand in wam_1b_6dof_1

% And it creates what's called a PanelGeo(metry), and this is what actually
% creates the .gdf file. We don't have to point the FloatingBody to a .gdf
% file like we did in wam_1b_6dof_1.

% We can look at the Geometry..
figure;
subplot(1,2,1);
plot(cyl.PanelGeo);
axis equal;
title({'Cylinder Panelization', 'dia: 10 m, draft: 10 m, height: 15 m',...
    'Body Coordinates'});

% THE VISUALISATION DOESN'T SEEM TO SHOW WHAT THE BELOW COMMENTS SUGGEST,
% BUT THE GDF FILE LOOKS TO BE DEFINED CORRECTLY NONETHELESS.
% The result may look a little strange at first.. The bottom is at z = -7.5
% and the top is at z = 2.5. So what's happening. 1) Only the part below
% the free surface is panelized. Here, the draft is 10, so we have panels
% of length 10 in the z direction (2.5 - -7.5). 2) It's drawn in body
% coordinates, the body origin is at (0,0,0).  When it goes into Wamit, we
% tell Wamit the Z position of the origin and it effectivly shifts it by
% its z position. 

globCylPanGeo = PanelGeo(cyl.PanelGeo);
globCylPanGeo.Translate([0 0 cyl.Zpos]);
subplot(1,2,2);
plot(globCylPanGeo);
axis equal;
title('Global Coordinates');

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
wam_run.Beta = 0;       % set the incident wave direction(s) (in radians)
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

wam_run.IReadRAO = 1; % 0 - use .4 file RAOs, 1 - use user-input RAOs.

wam_run.FloatingBodies = cyl;       % Now set the floating body we just 
                                    % created as the body used in Wamit run
wam_run.WriteRun;                   % this writes all the necessary Wamit 
                                    % input files
                                    
%% Run Wamit

% The following are default values, but just to show that they can be set
wam_run.ExePath = 'N:\wamitv7';           % This points to the location 
                                            % of the Wamit.exe
wam_run.ScratchPath = 'N:\wamitv7\scratch'; % Wamit needs a scratch folder
wam_run.UseridPath = 'N:\wamitv7';          % Location of UserId (license)

wam_run.Run;                           % Runs Wamit on the run that   
                                            % just created. Ties up Matlab
                                            
% wam_run.Run('Background');             % Runs Wamit in a command
                                            % window. Matlab can be used,
                                            % but make sure you wait to
                                            % read the results until the
                                            % run is finished.

%% Read results

% the WamitResult class reads in Wamit output files. 
wam_result = WamitResult(wam_run);  % create with the WamitRunCondition

% wam_result = WamitResult();         % Or create by setting the folder and
% wam_result.Folder = folder;         % run name
% wam_result.RunName = run_name;

wam_result.ReadResult;              % Read the Wamit output files

%% Compute RAOs externally to WAMIT
% 1) Get the frequency domain forces computed by WAMIT.
% i.e. the added mass, hydrodynamic damping, excitation force, and
% hydrostatic stiffness matrix
freqDomForces = wam_result.FreqDomForces;

% 2) Create the constraint matrix to constrain the bodies together.
P_6dof = ConstraintMatComp.FixedBodies([0 0 0]); % This should simply give the identity matrix in the case of a single body
P = P_6dof(logical(cyl.Modes.Vector),...
logical(cyl.Modes.Vector)); % Remove rows and columns for DoFs not studied.

% 3) Create the Frequency Domain computation - type of IEnergyComp
freqDomComp = FreqDomComp(freqDomForces, [cyl], 'Constrained', P);

% 4) Add desired changes to damping and stiffness matrices
Dpto_6dof = zeros(6,6); Dpto_6dof(3,3) = 10000; 
Dpto = Dpto_6dof(logical(cyl.Modes.Vector),...
    logical(cyl.Modes.Vector));
freqDomComp.SetDpto(Dpto);

% 5) Derive motion RAOs
motionRAOs = freqDomComp.Motions;

%% Write RAOs to .rao file for input to WAMIT and alter .CFG file
wam_run.MotionRAOs = motionRAOs;
wam_run.WriteRunIReadRAO;

%% Rerun WAMIT to obtain drift forces
wam_run.Run;

%% Read results and plot drift forces
wam_result = WamitResult(wam_run);  % create the updated WamitRunCondition
wam_result.ReadResult; % Read the updated WAMIT result

driftForces = wam_result.DriftForces; % Extract drift forces

% Plot drift forces for all 6 DoFs of the cylinder and for all methods
saveLoc = 'C:\Users\AlfredCotten\Desktop\Mean drift force WAMIT trials\cylinder';
saveNames = {'_csf','_mom','_prs'};
legendNames = {'csf','mom','prs'};
dofsTitle = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};

% Plot different methods together on same figures
figure
for i = 1:6
    subplot(3,2,i)
    hold on
    for j = 1:size(driftForces,2)-1
        plot(2*pi./driftForces(j).T(:),abs(driftForces(j).meanDriftForces(:,i)),'-')
    end
    xlabel('Wave frequency / rad/s')
    ylabel('Force/Moment')
    title(dofsTitle{i})
end
legend(legendNames)
savefig([saveLoc '\driftForces_allMethods'])

% Plot different methods separately on different figures
for j = 1:size(driftForces,2)
    figure
    for i = 1:6
        subplot(3,2,i)
        plot(driftForces(j).T,abs(driftForces(j).meanDriftForces(:,i)))
        xlabel('Wave period / s')
        ylabel('Force/Moment')
        title(dofsTitle{i})
    end
    sgtitle(driftForces(j).methodName)
%     savefig([saveLoc '\driftForces' char(saveNames(j))])
end
    
% the FreqDomForces class holds information on: 
hydroForces = wam_result.FreqDomForces;
A = hydroForces.A;      % added mass (Nperiod x DoF x Dof), real value
B = hydroForces.B;      % hydrodynamic damping (Nperiod x DoF x Dof), real
C = hydroForces.C;      % hydrostatic restoring force (DoF x DoF), real
Fex = hydroForces.Fex;  % excitation force (Nperiod x Ndir x DoF), complex

T = hydroForces.T;      % wave periods

figure;
subplot(3,1,1);
plot(T, squeeze(A(:,3,3)));
title({'Results for Cylinder (dia: 10 m, draft: 10 m, height: 15 m)', ...
    'Added Mass'});
ylabel('kg');
subplot(3,1,2);
plot(T, squeeze(B(:,3,3)));
title('Hydrodynamic Damping')
ylabel('Ns/m');
subplot(3,1,3);
plot(T, abs(squeeze(Fex(:,1,3))));
title('Excitation Force');
xlabel('Period (s)');
ylabel('N');


%% Analyze results

% the FreqDomComp class computes values such as motions and power. 
% It takes a FreqDomForces and a FloatingBody at its inputs
hydroComp = FreqDomComp(hydroForces, cyl);

% power computation requires a linear mechanical damping. 
dof = cyl.Modes.DoF;
Dpto = zeros(dof, dof);     % PTO damping matrix of size DoF x DoF
Dpto(3,3) = 10^3;           % Set a damping in the heave mode of motion

hydroComp.SetDpto(Dpto);    % Set the damping on the HydroComp

RAO = hydroComp.Motions;    % The motions for a unit amplitude wave
                            % (Nperiod x Ndir x DoF), complex valued

power = hydroComp.Power;    % Get the computed power. The output is a 
                            % matrix of size Nperiod x Ndir x Dof. In this
                            % case, because the PTO damping is zero for all
                            % values except (heave, heave), the output
                            % power is zero in all DoF except heave.
                            
power = squeeze(power(:,1,3));  % Just get the non-zero, heave power

% Maybe we want relative capture width
% The following function returns the wave energy flux for unit amplitude
% waves through a 1 meter section. 
incFlux = IWaves.UnitEnergyFlux(rho, T, hydroComp.H).';    

RCW = power./incFlux./diameter;

figure;
subplot(3,1,1);
axe = plotyy(T, abs([squeeze(RAO(:,1,1)) squeeze(RAO(:,1,3))]), T, ...
    180/pi*abs(squeeze(RAO(:,1,5))));
ylabel(axe(1), 'Surge, Heave (m/m)');
ylabel(axe(2), 'Pitch (deg/m)');
title({'Results for Cylinder (dia: 10 m, draft: 10 m, height: 15 m)', ...
    'Response Amplitude Operator'});
legend('Surge', 'Heave', 'Pitch', 'Location', 'NorthWest');

subplot(3,1,2);
plot(T, power./1000);
title('Power in 2 m high waves')
ylabel('kW');

subplot(3,1,3);
plot(T, RCW);
title('Relative capture width');
xlabel('Period (s)');
ylabel('RCW');

