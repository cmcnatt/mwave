% script to make an energy flux computation on a floating cylinder

% to run this script you require
%% 

% 1. WAMIT
% 2. mwave (https://github.com/cmcnatt/mwave)

% CELL ORDER
% 1. create and plot body
% 2. compute hydrodynmaics in WAMIT
% 3. compute energy flux
% 4. plot metrics
% 5. compare methods

clear; clc; close all;

%% make a cylinder using mwave FloatingCylinder class
% uses Wamit low order panel method

length = 2;
beam = 30;
draft = 24.9;


plotArgs = {'view', [-60 -30], ...
    'xlim', [-6 6], 'ylim', [-6 6], 'zlim', [-11 3]};

N = 2;  

Nx = N;          % number of panels in axial
Ny = N*12;           % number of panels in beam
Nz = N*6;          % number of panels in depth

rho = 1025;         % fluid density

% FloatingFlap creates panel geo then computes and sets the
% mass properties
flap_flat = FloatingFlap(rho, length, beam, draft, Nx, Ny, Nz);

figure;
subplot(1,2,1);
plot(flap_flat.PanelGeo); 
axis equal
grid on
box on;
set(gca, plotArgs{:});
title('Complete flap');

subplot(1,2,2);
plot(flap_flat.PanelGeo, 'onlywet', 'ShowNorm', 'noint'); 
axis equal
grid on
box on;
set(gca, plotArgs{:});

title('Flux surface showing normals');

%% Set up and run WAMIT, and read results

wamitPath = [energyFluxToolPath '\wamitRun'];   % Folder where WAMIT files go 
if 0 == exist(wamitPath, 'dir')
    mkdir(wamitPath);
end
delete([wamitPath '\*']);

wamRun = WamitRunCondition(wamitPath, 'test');   % make WAMIT run object, setting the run folder
wamRun.Rho = rho;                           % Set the fluid density
wamRun.T = 2:0.2:12;                        % wave periods to run
wamRun.Beta = 0;                            % incident wave direction
wamRun.H = 25;                              % water depth

% set the modes of motion: [1 0 1 0 1 0] means surge, heave, pitch, (but not sway, roll, yaw)
flap_flat.Modes = ModesOfMotion([0 0 0 0 1 0]);  
wamRun.FloatingBodies = flap_flat;                % Set the flapinder to be run by WAMIT

% Tell WAMIT to compute pressure and velocity at points on the GDF panels
wamRun.ComputeBodyPoints = true;
wamRun.ComputeVelocity = true;

% RESOURCE
wamRun.NCPU = 6;
wamRun.RAMGBmax = 8;

wamRun.WriteRun;                            % Write the input files for WAMIT
wamRun.Run;                                 % Use matlab to run WAMIT

wamResult = WamitResult(wamRun);            % Create a WAMIT result object
wamResult.ReadResult;                       % Read the WAMIT output files

%% Compute energy flux body surf

waveBody = wamResult.WaveBody;              % Get the BodySurfWaveField object
freqDomForces = wamResult.FreqDomForces;    % Get the frequency-domain forces computed by WAMIT  

comp = FreqDomComp(freqDomForces, flap_flat);     % Create a frequency domain computation
xi0 = comp.Motions;                         % Use the frequency-domain computation to compute the 
                                            % complex-valued motion amplitudes.
                                            % Array of size Nt x Nbeta x DoF where
                                            %   - Nt: number of periods
                                            %   - Nt: number of wave directions
                                            %   - DoF: degrees of freedom
                                            
waveBody.BodyMotions = xi0;                 % Set those amplitudes on the waveBody

% The BodySufWaveField object can be used to get pressure values over the
% body surface: 
% press = waveBody.Pressure('Total');

% Compute the energy flux over the body surface:
% The diffracted wave field i.e. the wave field over the body held fixed
fluxD = waveBody.EnergyFlux('Diffracted');

% The total wave field i.e. the moving body (total = sum of radiated and
% diffracted). This is computed the body without any mechanical
% constraints, e.g. PTO damping
flux0 = waveBody.EnergyFlux('Total');

% Set a PTO damping value on the frequency domain computation
Dpto = 11.6101^9;
comp.SetDpto(Dpto);
xiP = comp.Motions;                     % Compute the motions with the PTO damping
pow = comp.PowerRAO;                    % Compute the Power RAO in kW/m^2
[peakPow, iPeak] = max(pow);            % get the peak power and index

waveBody.BodyMotions = xiP;             % Set the motions on the waveBody

% Get the flux with PTO damping
fluxP = waveBody.EnergyFlux('Total');

% fluxD, flux0, and fluxP are cell arrays of PanelGeo objects, where the
% values are set with the real values of the energy flux at each panel

% Now compute the wave energy flux surface for a spectrum
Hs = 3.5;
Tp = 7.5;
spec = Bretschneider(Hs, Tp, wamRun.T); % Create a Bretschneider spectrum object
fluxS = waveBody.EnergyFlux('Total', 'spectra', spec);
fluxS = fluxS{1};

avgPow = comp.AveragePower(spec);       % Average power of the flapinder WEC 
                                        % in the wave spectrum (kW)

%% Plot the motion and power RAOs 
% (to give a better understanding of device behaviour) 

linCol = [0, 0.4470, 0.7410];

figure;
ymul = [1 1 180/pi];
ylab = {'Pitch [deg/m]'};
for n = 1
    subplot(2,1,n);
    plot(comp.T, ymul(n)*abs(xi0(:,1,n)), '--', 'color', linCol);
    hold on;
    plot(comp.T, ymul(n)*abs(xiP(:,1,n)), 'color', linCol);
    grid on
    xlabel('Period [s]');
    ylabel(ylab{n});
    if n == 1
        title('RAO''s');
        legend('No PTO damping', 'With PTO damping');
    end
end

% Plot the power RAO
subplot(2,1,2);
plot(comp.T, pow);
grid on;
xlabel('Period [s]');
ylabel('Power [kW/m^2]');
hold on;
plot(comp.T(iPeak), pow(iPeak), 'x', 'color', linCol);

%% Plot the energy flux on the geometries

clim = 2*[-1 1];

figure;
subplot(2,2,1);
surf(fluxD{iPeak}./1000);
axis equal
grid on;
box on;
set(gca, 'clim', clim, plotArgs{:});
title('Diffracted Energy Flux');
caxis = colorbar;
pos = get(caxis, 'position');
set(caxis, 'position', pos);
ylabel(caxis, 'Energy Flux Density [kW/m^2]');

subplot(2,2,2);
surf(flux0{iPeak}./1000);
axis equal
grid on;
box on;
set(gca, 'clim', clim, plotArgs{:});
title('Energy Flux - No Damping');

subplot(2,2,3);
surf(fluxP{iPeak}./1000);
axis equal
grid on;
box on;
set(gca, 'clim', clim, plotArgs{:});
title('Energy Flux - Power Absorption');

subplot(2,2,4);
surf(fluxS./1000);
axis equal
grid on;
box on;
set(gca, 'clim', clim, plotArgs{:});
title('Energy Flux - Spectrum');

%% Table comparing power generation for each method

% Integrate over body surface
efD = fluxD{iPeak}.Values;      % Energy flux is in the Values property on the PanelGeo
ef0 = flux0{iPeak}.Values;      
efP = fluxP{iPeak}.Values;      
efS = fluxS.Values;
areas = fluxD{iPeak}.Areas;     % Get the area of each panel - this is actually the same
                                % for any of the PanelGeos or any index
energyFluxD = sum(efD.*areas)./1000;    % Diffracted energy flux
energyFlux0 = sum(ef0.*areas)./1000;    % 0 damping energy flux
energyFluxP = sum(efP.*areas)./1000;    % Power absorping energy flux
energyFluxS = sum(efS.*areas)./1000;     % Spectral energy flux

expectedPow = round([0 0 peakPow avgPow]', 2);
energyFlux = round([energyFluxD, energyFlux0, energyFluxP, energyFluxS]', 2);

tab = table(expectedPow, energyFlux)
