% This test sets up and runs Wamit on a 6 dof cylinder. 
% The cylinder's geometry file is made in matlab with mwave. This example 
% is the same as example_Wamit_1bod_6dof_1, except for that aspect. 

clear; close all; clc;

%% Set up the run

% every run needs a name. All the Wamit input and output files are named 
% based off this name
run_name = 'wam_1b_6dof_2';         

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

% The FloatingCylinder class computes all of the mass properties that we
% had to compute by hand in wam_1b_6dof_1

% And it creates what's called a PanelGeo(metry), and this is what actually
% creates the .gdf file. We don't have to point the FloatingBody to a .gdf
% file like we did in wam_1b_6dof_1.

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

wam_run.FloatingBodies = cyl;       % Now set the floating body we just 
                                    % created as the body used in Wamit run
wam_run.WriteRun;                   % this writes all the necessary Wamit 
                                    % input files
                                    
%% Run Wamit

wam_run.Run;                            % Runs Wamit on the run that was 
                                        % just created. Ties up Matlab
                                            
% wam_run.Run('Background');            % Runs Wamit in a command window. 
                                        % Matlab can be used, but make sure 
                                        % you wait to read the results 
                                        % until the run is finished.
                                        
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
grid on;
subplot(3,1,2);
plot(T, squeeze(B(:,3,3)));
title('Hydrodynamic Damping')
ylabel('Ns/m');
grid on;
subplot(3,1,3);
plot(T, abs(squeeze(Fex(:,1,3))));
title('Excitation Force');
xlabel('Period (s)');
ylabel('N');
grid on;


%% Analyze results

% the FreqDomComp class computes values such as motions and power. 
% It takes a FreqDomForces and a FloatingBody at its inputs
hydroComp = FreqDomComp(hydroForces, cyl);

% power computation requires a linear mechanical damping. 
dof = cyl.Modes.DoF;
Dpto = zeros(dof, dof);     % PTO damping matrix of size DoF x DoF
Dpto(3,3) = 5*10^4;           % Set a damping in the heave mode of motion

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
grid on;

subplot(3,1,2);
plot(T, power./1000);
title('Power in 2 m high waves')
ylabel('kW');
grid on;

subplot(3,1,3);
plot(T, RCW);
title('Relative capture width');
xlabel('Period (s)');
ylabel('RCW');
grid on;

