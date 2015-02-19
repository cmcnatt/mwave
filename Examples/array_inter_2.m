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
% In this example, we'll use the HydroBody that we created in example
% Wamit_createHB_2 to compute power in an array from a spectral wave field.
%
% For more info on Array computations, see example array_inter_1, and for
% more info on HydroBodies see example hydroBody_1

%% Load the HydroBody

% load 'wam_hb2_1_hb';    % Load the HydroBody based on its name..

% Or we can be more explicit if necessary         
load([mwavePath 'Examples\HydroBodies\wam_hb2_1_hb']);

%% Create an the the initial array layout and array computation object
% The body positions and incident waves can be changed later

% 5 bodies at the following positions
%
%   *(0,10)
%
%   *(0,5)
%
%   *(0,0)
%
%   *(0,-5)
%
%   *(0,-10)

% Create an array of HydroBodies
for m = 1:5
    hbs(m) = HydroBody(hydBody);    
end

% The p
hbs(1).XYpos = [0, 10];
hbs(2).XYpos = [0, 5];
hbs(3).XYpos = [0, 0];
hbs(4).XYpos = [0, -5];
hbs(5).XYpos = [0, -10];

% Initialize the array computation with unit amplitude plane waves
T = hydBody.T;         
h = hydBody.H;      
a = ones(size(T)); 
beta = 0;

iwaves_ua = PlaneWaves(a, T, beta, h);

% create the initial array computation object
arrayComp = HydroArrayComp(hbs, iwaves_ua);

%% Let's take a quick look at the body and its frequency response: RAO

geo = hydBody.PanelGeo;

figure;
subplot(3,1,1);
plot(geo);
axis equal;
title('Single Cylinder');

hComp = HydroBodyComp(hydBody, iwaves_ua);

nondimDpto = 2.94;      % Value from Child and Venugopal (2010)
rho = hydBody.Rho;
a = 1;
a2k0 = 0.8;
omega0 = 2*a./a2k0;
dpto = nondimDpto*rho*a^3*omega0;

hComp.SetDpto(dpto);
hComp.SetK(0);          % You can also use this to set the reactive tuning
xi = hComp.Motions;

pow1 = hComp.Power;

subplot(3,1,2);
plot(T, abs(xi));
xlabel('T (s)');
ylabel('RAO (m/m)');

subplot(3,1,3);
plot(T, 1/1000*pow1);
xlabel('T (s)');
ylabel('Power (kW)');

%% HydroBody Computations for directional spectrum
% This example uses a fully directional spectrum. The other two examples
% will use more simplified representations of a wave spectrum: a
% unidirectional spectrum and a single plance wave. It is interesting to
% compare the results of the power absoption

% create a wave spectrum
Hs = 2;                         % significant wave height
Tm = 4;                         % modal period 
f = 1./T;                       % use the wave periods defined by the 
                                % HydroBody. If the periods are ascending,
                                % then the frequencies will be decending,
                                % but that's ok. The order must be the same
                                % as it is in the HydroBody
S = bretschneider(Hs, 1./Tm, f);% the bretschneider function is built into 
                                % mwave. For another spectrum (e.g.
                                % JONSWAP), you would need to come up the
                                % the spectram points yourself.  

% next create a directional spread
s = 6;                              
betac = 0;                          
beta = linspace(-pi/2, pi/2, 9);    
G = cosSpectSpread(s, betac, beta); 
SG = S'*G;                          

% This is the wave spectrum that will be used to define the incident
% amplitudes when computing power in the array
Spec = WaveSpectrum(SG, f, beta);   

% Let's take a quick look at the omni directional spectrum and RAO on the
% same axis as well as the directional spectrum;
figure;
subplot(2,1,1);
ax = plotyy(T, Spec.Spectrum('Nondir'), T, 1/1000*pow1);
ylabel(ax(1), 'Spectral density (m^2/Hz)');
ylabel(ax(2), 'WEC power response (kW)');  % Note that the power response 
                                           % here is for unite amplitude 
                                           % waves, not for the waves
                                           % defined by the spectrum!

subplot(2,1,2);
S2 = Spec.Spectrum;
surf(T, beta, S2');
xlabel(gca, 'T (s)');
ylabel(gca, 'Dir (radian)');
zlabel(gca, 'Spectral density (m^2/Hz/rad)');
title(gca, 'Directional spectrum');


a = Spec.Amplitudes;    % Get the amplitude of the incident waves from the 
                        % wave spectrum. This is how we get the wave
                        % spectrum into the array compuation

for m = 1:length(beta)
    iwaves_ds(m) = PlaneWaves(abs(a(:,m)), T, beta(m), h);
end

% Now we set the incident waves of the array computation, which came frome
% the wave spectrum
arrayComp.IncWaves = iwaves_ds;

%% Compute power for array in spectrum
% set the PTO damping 
% dpto comes from the section above
dof = arrayComp.DoF;    % Dof for entire array
Dpto = zeros(dof,dof);
for m = 1:5     
    Dpto(m,m) = dpto;   % in this case, because the cylinder only has one 
                        % DoF, the damping values are the entire diagonal
end

arrayComp.SetDpto(Dpto);    

% The power is a matrix of size: nT x nBeta x nBody
power = arrayComp.Power; 

% total power for the entire spectrum is the sum of all the power values
% from each wave component on each body
Ptot_dirSpec = sum(sum(sum(power)));

% the power absorbed per body
Pbod = squeeze(sum(sum(power)));

% the total power absorbed for all bodies as a function of frequency
Pfreq_dirSpec = squeeze(sum(squeeze(sum(power, 2)), 2));

figure;
plot(T, 1/1000*Pfreq_dirSpec);
xlabel('T (s)');
ylabel('Power absorbed (kW)');
title('Power for 5 body array in directional seas'); 

%% HydroBody Computations for unidirectional spectrum

% create a wave spectrum - this is the same spectrum as above, but without
% directional spreading
Hs = 2;                         
fm = 0.2;                       
f = 1./T;                       
S = bretschneider(Hs, fm, f);   

Spec = WaveSpectrum(S, f);   

a = Spec.Amplitudes;    
beta = 0;

iwaves_us = PlaneWaves(abs(a), T, beta, h);

% Now we set the incident waves of the array computation, which came frome
% the wave spectrum
arrayComp.IncWaves = iwaves_us;

% PTO damping is already set, but can be changed...
power = arrayComp.Power; 

% Total power for array
Ptot_uniSpec = sum(sum(power));

% the total power absorbed for all bodies as a function of frequency
Pfreq_uniSpec = squeeze(sum(squeeze(sum(power, 2)), 2));

%% HydroBody Computations for 'equivalent wave'

% define a wave spectrum by a single plane wave that is 'equivalent' to the
% wave spectrum
Hs = 2;                         
fm = 0.2; 

a = zeros(size(T));                             % Amps for all wave periods 
[~, ifm] = min(abs(T - (1/fm)*ones(size(T))));  % are zero except at that 
a(ifm) = Hs/2;                                  % of the period closest to
                                                % the modal period
beta = 0;

iwaves_ew = PlaneWaves(a, T, beta, h);

arrayComp.IncWaves = iwaves_ew;

power = arrayComp.Power; 

% Total power for array
Ptot_eWave = sum(sum(power));

%% Compare the total power from the three different wave representations

[Ptot_dirSpec, Ptot_uniSpec, Ptot_eWave]

%% Go back to the directional spectrum incident waves and change the body 
% positions and recompute

arrayComp.IncWaves = iwaves_ds;

%  5 bodies, (x, y) positions
pos = zeros(5, 2);

%               *(5,10)
%
%   *(-5,5)
%
%               *(5,0)
%
%   *(-5,-5)
%
%               *(5,-10)

pos(1,:) = [5, 10];
pos(2,:) = [-5, 5];
pos(3,:) = [5, 0];
pos(4,:) = [-5, -5];
pos(5,:) = [5, -10];

arrayComp.BodXY = pos;

power = arrayComp.Power;

Ptot_newPos = sum(sum(sum(power)));

% compare with original position

[Ptot_dirSpec, Ptot_newPos]

%% Take a look at the wave field by looking at the significant wave height

isarray = true;
x = -20:0.5:20;
[X, Y] = meshgrid(x, x);    % square

waveField = arrayComp.WaveField(isarray, X, Y, 'NoVel'); 

Hs_wf = waveField.SigWaveHeight('Total');

pcolor(X,Y,Hs_wf./Hs);
fet;
set(gca, 'clim', [0.95 1.05]);
colorbar;
title('Disturbance coefficient (Hs/Hs0)');