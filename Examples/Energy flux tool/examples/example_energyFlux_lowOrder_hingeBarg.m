% script to make an energy flux computation on a floating cylinder

% to run this script you require

% 1. WAMIT
% 2. mwave (https://github.com/cmcnatt/mwave)

% CELL ORDER
% 0. Set up
%   Select one of 1a, 1b
% 1. import low order GDF
% 2. compute hydrodynmaics in WAMIT
% 3. compute energy flux
% 4. plot metrics
% 5. compare methods

clear; clc; close all;

%% 0. General set up

wamitPath = [energyFluxToolPath '\wamitRun'];   % Folder where WAMIT files go 
if 0 == exist(wamitPath, 'dir')
    mkdir(wamitPath);
end
delete([wamitPath '\*']);

% device = '\moceanGdfs';
% gdfPath = [energyFluxToolPath '\geoFiles' device];     % Folder of geometry files
rho = 1025;                                     % fluid density

plotArgs = {'view', [-40 30], ...
    'xlim', [-30 30], 'ylim', [-6 6], 'zlim', [-5 3]};

run_name = 'barge';  

%% 1. Make a Barge Raft geometry using MWave

rho = 1025;
% the barge hull parameters
len = 30;
beam = 8;
draft = 2;
hei = 4;
space = 1;

hingeLoc = [0, 0, 0];

% make barge

Ny = beam; Nx = len; Nz = draft;

cg_fwd = [-(len + space)/2, 0, 0];
fwd = FloatingBox(rho, len, beam, hei, draft, Nx, Ny, Nz,'Translate', cg_fwd);

% Create the floating bodies
% wamit panel size
wamPanSize = 4;

fwd.Cg = [-(len + space)/2, 0, 0];                                    
fwd.Modes = ModesOfMotion([1 0 1 0 1 0]);  
fwd.ISurfPan = 0;
fwd.WamILowHi = 0;
fwd.SurfAboveZ0 = true;
fwd.WamPanelSize = wamPanSize;


cg_aft = [(len + space)/2, 0, 0];
aft = FloatingBox(rho, len, beam, hei, draft, Nx, Ny, Nz,'Translate', cg_aft);
aft.Cg = [(len + space)/2, 0, 0];                                    
aft.Modes = ModesOfMotion([1 0 1 0 1 0]);  
aft.ISurfPan = 0;
aft.WamILowHi = 0;
aft.SurfAboveZ0 = true;
aft.WamPanelSize = wamPanSize;

% plot the geometry

figure;
surf(fwd.PanelGeo, 'Color', MColor.Yellow);
hold on;
surf(aft.PanelGeo, 'Color', MColor.Yellow);
axis equal
grid on
set(gca, 'view', [20 10]);
title('Hinged Barge Raft');


%% 2. Set up and run WAMIT, and read results

% Set up and run WAMIT
wamRun = WamitRunCondition(wamitPath, run_name); 
wamRun.UseDirectSolver = 1;
wamRun.NCPU = 6;
wamRun.RAMGBmax = 15;

wamRun.Rho = rho;      
wamRun.T = 4:0.2:12;   
%wamRun.T = 4:0.5:8;   
wamRun.Beta = 0;       
wamRun.H = Inf;           

% Tell WAMIT to compute pressure and velocity at points on the GDF panels
wamRun.ComputeBodyPoints = true;
wamRun.ComputeVelocity = true;

wamRun.FloatingBodies = [fwd, aft];       
wamRun.WriteRun;                  

wamRun.Run;   

wamResult = WamitResult(wamRun);  
wamResult.ReadResult; 


%% 3. Compute energy flux body surf

waveBody = wamResult.WaveBody;              % Get the BodySurfWaveField object
freqDomForces = wamResult.FreqDomForces;    % Get the frequency-domain forces computed by WAMIT  

% Create the constraint matrix to constrain the bodies together.
origin = [0 0 0];                         % Origin about which you want 
                                          % the motions to be referenced
P = ConstraintMatComp.HingedBodies([fwd.Cg; aft.Cg], hingeLoc, 'Origin', origin, 'Planar');

comp = FreqDomComp(freqDomForces, [fwd, aft], 'Constrained', P);     % Create a frequency domain computation
%xi0 = comp.Motions;                         % Use the frequency-domain computation to compute the 
xi0=comp.Motions('orgcoor');
                                            % complex-valued motion amplitudes.
                                            % Array of size Nt x Nbeta x DoF where
                                            %   - Nt: number of periods
                                            %   - Nbeta: number of wave directions
                                            %   - DoF: degrees of freedom
   
waveBody.BodyMotions = xi0;                 % Set those amplitudes on the waveBody

xi0=comp.Motions('orgcoor');

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

%% Set a PTO damping value on the frequency domain computation
Dpto = zeros(4,4);                          
Dpto(4,4) = 10^6;
comp.SetDpto(Dpto);
xiP = comp.Motions('orgcoor');          % Compute the motions with the PTO damping
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

avgPow = comp.AveragePower(spec)       % Average power of the cylinder WEC 
                                        % in the wave spectrum (kW)

%% 4. Plot the motion, power RAOs, and flux maps 
% (to give a better understanding of device behaviour) 

linCol = [0, 0.4470, 0.7410];

figure;
ymul = [1 1 180/pi];
ylab = {'Surge [m/m]', 'Heave [m/m]', 'Pitch [deg/m]'};
for n = 1:3
    subplot(4,1,n);
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
subplot(4,1,4);
plot(comp.T, pow);
grid on;
xlabel('Period [s]');
ylabel('Power [kW/m^2]');
hold on;
plot(comp.T(iPeak), pow(iPeak), 'x', 'color', linCol);

%% Plot the energy flux on the geometries
clim = 5*[-1 1];

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

%% 5. Table comparing power generation for each method

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