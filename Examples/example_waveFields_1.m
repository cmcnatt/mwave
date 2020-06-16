% This example shows some of the basics of the WaveField and
% WaveFieldCollection classes

clear; close all; clc;

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

%% WaveFieldCollection example
% A WaveFieldCollection is a way to store multiple related wave fields in
% a single object. A common way to use a WaveFieldCollection is to store 
% multiple wave directions at the same period. A single WaveField does not
% allow one to have multiple directions at the same period, because these
% waves would be coherent. However, in a spectral sense, these wave
% components are orthongal and so it is justifiable to keep them separate.
%
% In this example a wave field will be generated with a directional
% spectrum

disp('ERROR: This example section is not working correctly')

% create a wave spectrum
% first create a frequency point spectrum
Hs = 2;                         % significant wave height
fm = 0.2;                       % modal period
f = 0.08:0.02:0.8;              % the range of wave frequencies (Hz)
S = bretschneider(Hs, fm, f);   % Bretschneider spectrum

% next create a directional spread
s = 4;                              % spreading parameter
betac = 0;                          % center direction
beta = linspace(-pi, pi, 21);       % directions (radians)
G = cosSpectSpread(s, betac, beta); % spreading function
SG = S'*G;                          % create a matrix of the spectral amps

% finally, create the wave spectrum object
Spec = WaveSpectrum(SG, f, beta);   

a = Spec.Amplitudes('RandPhase');   % return wave amplitudes with random 
                                    % phase for each wave component in the
                                    % spectrum
                                    
f = Spec.Frequencies;               % Already had f,beta from above, but 
beta = Spec.Directions;             % just wanted to show these methods
S = Spec.Spectrum;                  % These are the actual Spectral values 
                                    % (it is the same as the SG variable
                                    % above that was used to create the
                                    % WaveSpectrum object)
                   
figure;                             % plot the spectrum
surf(f, 180/pi*beta, S.');
xlabel('frequency (Hz)');
ylabel('direction (deg)');
zlabel('Spectral density (m^2/Hz)');

% create wave field for this spectrum. All the waves at the same direction 
% can be grouped into a single WaveField. Then the WaveFields at each
% direction are grouped together in a WaveFieldCollection
rho = 1025;
h = 40;
x = -50:0.5:50;
y = -50:0.5:50;

[X, Y] = meshgrid(x, y);

for n = 1:length(beta)
    wcs = PlaneWaves(abs(a(:,n)), 1./f, beta(n)*ones(size(f)), h, ...
        angle(a(:,n)));
    wfs(n) = PlaneWaveField(rho, wcs, 1, X, Y);
end

% creates the WaveFieldCollection object
wf = WaveFieldCollection(wfs, 'Direction', beta);   % a WaveFieldCollection 
                                                    % has what is called a
                                                    % collection type. This
                                                    % is specified by the
                                                    % second argument:
                                                    % 'Direction', and it
                                                    % has what is called an
                                                    % index, which is
                                                    % specified by the last
                                                    % argument, beta

wf.CollType     % The CollType should be 'Direction'
wf.WFcount      % The WFcount should be the same as the number of direction
wf.Indices      % The Indices are the directions

% A WaveFieldCollection is a WaveField and can output regular things like
% wave elevation, which is a cell array of plane wave field
eta = wf.Elevation; 

figure;             % lets just look at a few that have wave energy
istartF = 5;
istopF = 9;
nf = istopF - istartF + 1;
istartB = 9;
istopB = 13;
nb = istopB - istartB + 1;

for m = istartB:istopB
    for n = istartF:istopF
        subplot(nf, nb, (m - istartB)*nb + n - istartF + 1);
        pcolor(X,Y,real(eta{n, m}));
        if(m == istartB)
            title(['f = ' num2str(f(n))]);
        end
        if (n == istartF)
            ylabel(['Dir = ' num2str(beta(m))]);
        end
        set(gca, 'xtick', [], 'ytick', [], 'clim', 0.1*[-1 1]);
        shading flat;
        axis equal;
        axis tight
    end
end

% SigWaveHeight - gives a wave field picture of the significat wave height
% In this case it is not very interesting, because Hs is the same
% everywhere. (It is more interesting when there are bodies present that
% modify the spectral wave field).
Hs_wf = wf.SigWaveHeight('Merge');
Hs_spec = Spec.SigWaveHeight;       % Hs can also be computed from the 
                                    % spectrum
                                    
% lets see how these values compare the Hs originally specified
[Hs Hs_wf(1,1) Hs_spec]

% there is some error due to how the spectrum is discritized into wave
% components, but it's pretty close
