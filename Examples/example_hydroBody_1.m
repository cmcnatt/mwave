% In this example, we'll use the HydroBody we created in example
% example_Wamit_createHB_1 to compute power from diffrent incident wave 
% direction and to compute wave fields.

%% Load the HydroBody

% load 'wam_hb1_1_hb';    % Load the HydroBody based on its name..

run_name = 'wam_hb_1';  % Or we can be more explicit if necessary         
folder = [mwavePath 'Examples\HydroBodies'];  
load([folder '\wam_hb1_1_hb']);

%% HydroBody Computations

% A HydroBody is a body (a FloatingBody) with hydrodynamic wave information
% on top of it. The wave info is:
%
%   - added mass
%   - hydrodynamic damping
%   - hydrostatic matrix
%   - force transfer matrix (like excitation force)
%   - radiation coefficients (creates the radiated waves)
%   - diffraction transfer matrix (creates the scattered waves)
%
% It is all dependent of wave period and water depth. A single HydroBody
% can hold info on multiple periods, but only a single water depth. For
% different depths, new HydroBodies need to be computed. (Again, Inf is not
% possible at this time). 
% 
% HydroBodies don't really do anything by themselves - they need incident
% waves to react to. These waves can come from any direction. 
%
% Computations are done with the handy FreqDomComp class, which has been
% used in other examples, except here it takes a different set of arguments
% in the constructor - a HydroBody and IWaves (The "I" in IWaves means 
% interface, not incident). We've loaded the HydroBody, now let's create a 
% set of incident plane waves.

T = hydBody.T;      % The waves have to have the same periods at the 
                    % HydroBody
h = hydBody.H;      % And the same water depth
beta = [0 pi/4 pi/3, pi/2]; % 4 incident wave directions
a = ones(size(T)); % Just unit amplitude waves, but could come from a 
                    % spectrum (or anything)

for n = 1:length(beta)
    % Here, we're building an array of plane waves at differnt directions.
    % Each set of plane waves can only have a single incident direction.
    % For more info see the example: waveFields_1
    iwaves(n) = PlaneWaves(a, T, beta(n), h);
end

% Now create our comptuation with the HydroBody and the PlaneWaves
hydroComp = FreqDomComp(hydBody, iwaves);

% We can use the HydroComp just like we did before to compute power, etc.
dof = hydBody.DoF;
Dpto = zeros(dof, dof);     
Dpto(7,7) = 10^8;           
Dpto(8,8) = 10^8;           
hydroComp.SetDpto(Dpto);    

power = hydroComp.Power;    
power = squeeze(power(:,:,7)) + squeeze(power(:,:,8));   % Here we have 4 
                                                         % different inc
                                                         % directions

figure;
plot(T, [squeeze(power(:,1)), squeeze(power(:,2)), squeeze(power(:,3)), ...
    squeeze(power(:,4))]./1000);
title('Power in 2 m high waves')
ylabel('kW');
xlabel('Period (s)');
Betas = {'\beta = 0', '\beta = \pi/4', '\beta = \pi/3', '\beta = \pi/2'};
legend(Betas);

% Now you might be thinking, what's the big deal, we could have just used
% Wamit to compute at different incident wave directions. And of course we
% have info on those incident directions, we computed the HydroBody with 40
% incident directions, remember. Well that is kind of the point, the
% excitation force information from those 40 directions is stored in the
% more compact "force transfer matrix". And the HydroBody also can compute
% WaveFields at arbitrary points! and can be used in array computations:
% see example array_inter_1

% So, let's see the wave fields. 
isarray = true;
x = -100:1:100;
[X, Y] = meshgrid(x, x);    % square

waveField = hydroComp.WaveField(isarray, X, Y, 'NoVel');    % creates an 
                                                            % FBWaveField.
                                                            % 'NoVel' means
                                                            % no velocity,
                                                            % which saves
                                                            % time

eta = waveField.Elevation('Total');             % already has motions

sect = hydBody.WaterPlaneSec;
thet = (0:(2*pi/50):2*pi)';
cir = hydBody.Rcir*[cos(thet) sin(thet)];
iT = 4;
figure;
for n = 1:4
    subplot(2,2,n);
    pcolor(X,Y,abs(eta{iT,n}));
    hold on;
    plot(sect(:,1), sect(:,2), 'w');
    plot(cir(:,1), cir(:,2));
    fet;
    set(gca, 'clim', [0.7 1.3]);
    title(Betas{n});
end

% The white outline is the outline of the body and the black circle is the
% circumscribing circle. Note how the wave field is just a big red then
% blue splotch near the center of the circle. This is because the theory
% does not produce accurate results inside the cylinder (although, the wave
% field actually is pretty good just inside the cylinder). 




