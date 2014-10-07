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
% This example shows some of the basics of the WaveField and
% WaveFieldCollection classes

%% Create a basic WaveField
% The most basic wavefield really doesn't compute much (well, it does
% compute evelvation, significant wave height, and energy flux). It is
% mainly a holder for pressure and velocity values that were computed by
% something else (e.g. Wamit).
%
% Here we'll make a "meaningless" wave field - i.e. the pressure and
% velocity values will be arbitrary values

rho = 1000;                 % Fluid density (kg/m^3)
g = IWaves.G;               % Gravitational constant - convenient to get it 
                            % from IWaves.G
h = 100;                    % Water depth (m)
t = 2:6;                    % Wave periods (s)

% first we'll create an array type WaveField
isarray = true;
x = -50:1:50;
y = -100:1:100;
[X, Y] = meshgrid(x, y);    % create the meshgrid of (x,y) spatial points

nT = length(t);
[nX, nY] = size(X);

% Pressure and velocity are complex! to represent the harmonic nature of 
% the linear wave field.
p = (1 + 1i)*ones(nT, nX, nY);      % the pressure values are in an array.                                    
vel = (1 + 1i)*ones(nT, 3, nX, nY); % the velocity values also in an array.
                                    % The 3 represents the (x,y,z) 
                                    % component of velocity

% and here is the wave field with zero values..
wf = WaveField(rho, g, h, t, p, vel, isarray, X, Y);

% For more interesting stuff, see the rest of this example..

% now, we'll make a point list type wave field.
isarray = false;
nP = 64;
theta = (0:2*pi/nP:2*pi*(1-1/nP))';             % the points will be on a 
points = [sin(theta), cos(theta), zeros(nP, 1)];% circle or radius 1 in 
                                                % (x,y) at z = 0;
                                                
p = (1 + 1i)*ones(nT, nP);          % an nT x nP array of complex values
vel = (1 + 1i)*ones(nT, 3, nP);     % an nT x 3 x nP array of complex

wf = WaveField(rho, g, h, t, p, vel, isarray, points);

% Now, what can we do with a wave field.. 
pOut = wf.Pressure;                 % get back our pressure values
velOut = wf.Velocity;               % get back our velocity
eta = wf.Elevation;                 % get elevation, which is p/(rho*g)

% These are the same values that we put in, in our arrays p, vel, but they
% are not in a Cell array (which I found convienent for handling the data).

%% Create a WaveField of plane waves 
% The last example created the most basic wave field where the pressure and
% velocity had been computed by another program (e.g. Wamit), but mwave
% also has some classes that can do some work. Classes build on
% FuncWaveField (Functional WaveField) can actually compute values of the
% wave field in space (it's just linear wave theory after all). 
%
% This example shows how to create a wave field of plane waves.

N = 12;                         % Number of wave components

rho = 1000;                     % Fluid density (kg/m^3)
a = ones(1,N);                  % The amplitude of each component (m)
t = linspace(1,12,N);           % The wave period (s) of each component
beta = 0;                       % The incident wave direction (radians!) 
                                % for the plane waves. All waves must have 
                                % the same direction in a PlaneWaveField 
                                % object. For multiple directions, use the
                                % WaveFieldCollection class - more on this
                                % below.
h = 40;                         % Water depth (m)
epsilon = 2*pi*rand(size(a));   % Phase

waves = PlaneWaves(a, t, beta, h, epsilon); % The PlaneWaves class is NOT
                                            % a WaveField. These are wave
                                            % components - they hold
                                            % infomation need to create the
                                            % WaveField
lambda = waves.Lambda;          % for example - the wave lengths

isarray = true;
x = -50:0.5:50;
[X, Y] = meshgrid(x, x);        % a square    

wf = PlaneWaveField(rho, waves, isarray, X, Y);

eta = wf.Elevation;

% Now let's take a look...
figure;
for n = 1:N
    subplot(3,4,n);
    pcolor(X,Y,real(eta{n}));   % Thre real() is the wave at time: t = 0
    fet;    % This is a shortcut for: shading flat; axis equal; axis tight
    title(['T = ' num2str(t(n)) ' s']);
    if (n == 9)
        xlabel('x (m)');
        ylabel('y (m)');
    else
        set(gca, 'xtick', [], 'ytick', []);
    end
    set(gca, 'clim', [-1 1]);
end

%% WaveField math
% Math can be performed on WaveFields with overloaded. Here we'll make a
% wave field that would be created under point absorber theory by adding
% together a PlaneWaveField (the incident wave) and a PntSrcWaveField (the
% radiated wave of a surging (in this case) point source). We can change
% the amplitude and phase of the point source motion and see how it affects
% the wave field

rho = 1000;                     
a = 1;
t = 6;
beta = 0;
h = 40;                         

waves = PlaneWaves(a, t, beta, h); 

x = -100:1:100;
[X, Y] = meshgrid(x, x);        

% This plane wave field is our incident wave field
iwf = PlaneWaveField(rho, waves, true, X, Y);

type = 'Surge';     % Surging point source (could be 'Heave' or 'Sway')
loc = [0 0];        % (x,y) location of the point source - here the origin

% This is our radiated wave field.
rwf = PntSrcWaveField('Surge', loc, rho, waves, true, X, Y);

% Let's see what the unadjusted WaveFields look like

etaI = iwf.Elevation;
etaR = rwf.Elevation;

figure
subplot(3,3,1);
pcolor(X, Y, real(etaI{1}));
title('Incident wave');
set(gca, 'clim', [-1 1], 'xtick', [], 'ytick', []);
fet;

subplot(3,3,2);
pcolor(X, Y, real(etaR{1}));
title('Unadjusted radiated wave');
set(gca, 'clim', [-1 1], 'xtick', [], 'ytick', []);
fet;

% OK, now for the math..

A = 1 + 1i;         % This is the complex amplitude to modify the radiated 
                    % wave field - it is like changing the amplitude and
                    % phase of a point source type WEC

rwf2 = A*rwf;       % (constant)*(WaveField)
twf = iwf + rwf2;   % (WaveField) + (WaveField); This could also be written
                    % as: twf = iwf + A*rwf

etaI = iwf.Elevation;
etaR = rwf2.Elevation;
etaT = twf.Elevation;

subplot(3,3,4);
pcolor(X, Y, real(etaI{1}));
title('Incident wave');
set(gca, 'clim', [-1 1], 'xtick', [], 'ytick', []);
fet;

subplot(3,3,5);
pcolor(X, Y, real(etaR{1}));
title('Adjusted radiated wave');
set(gca, 'clim', [-1 1], 'xtick', [], 'ytick', []);
fet;

subplot(3,3,6);
pcolor(X, Y, real(etaT{1}));
title('Total wave (Real)');
set(gca, 'clim', [-1 1], 'xtick', [], 'ytick', []);
fet;

subplot(3,3,9);
pcolor(X, Y, abs(etaT{1}));
title('Total wave (Abs)');
set(gca, 'clim', [0.7 1.3], 'xtick', [], 'ytick', []);
fet;

% Note that you can see the wave shadow behind the point source in the
% total wave field. You can adjust the 'A' value and see how it affects the
% wave field.

%% TODO: WaveFieldCollection example
