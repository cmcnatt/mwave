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
% This test sets up and runs Wamit on a 6 dof cylinder. 
% The cylinder's geometry file is made in matlab with mwave. This example 
% is the same as Wamit_1bod_6dof_1, except for that aspect. 

%% Set up the run

% every run needs a name. All the Wamit input and output files are named 
% based off this name
run_name = 'wam_1b_6dof_2';         

% This is the folder where all the Wamit input and output files go (You 
% need to make this folder!). It is good practise to name the folder the 
% same as the run name.
folder = [mwavePath '\Examples\BemRuns\' run_name];  
                                                    
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
cyl.Modes = ModesOfMotion([1 1 1 1 1 1]);   % This sets up the modes of 
                                            % motion for the cylinder.
                                            
% Create a WamitRunCondition object - this builds the Wamit run
wam_run = WamitRunCondition(folder, run_name);  

wam_run.Rho = rho;      % set the fluid density (kg/m^3)
wam_run.T = 4:0.1:12;   % set the wave period(s)(in s)
wam_run.Beta = 0;       % set the incident wave direction(s) (in radians)
wam_run.H = Inf;        % set the water depth (in m) - can be a positive 
                        % value or Inf     

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

% the FreqDomForces class holds information on: 
hydroForces = wam_result.FreqDomForces;
A = hydroForces.A;      % added mass (Nperiod x DoF xDof), real value
B = hydroForces.B;      % hydrodynamic damping (Nperiod x DoF xDof), real
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

